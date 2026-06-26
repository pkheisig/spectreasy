import { useEffect, useMemo, useRef, useState } from 'react';
import axios from 'axios';
import { Download, FileText, Plus, Upload } from 'lucide-react';
import './PanelBuilder.css';

const API_BASE = (() => {
    const envBase = (import.meta.env.VITE_API_BASE as string | undefined)?.trim();
    if (envBase) return envBase.replace(/\/$/, '');
    if (typeof window !== 'undefined') {
        if (window.location.port === '5174') return 'http://localhost:8000';
        return window.location.origin.replace(/\/$/, '');
    }
    return 'http://localhost:8000';
})();

type LibraryInfo = {
    id: string;
    label: string;
};

type DetectorInfo = {
    detector: string;
    label: string;
    laser: string;
    emission: number;
    color: string;
};

type FluorInfo = {
    fluorophore: string;
    peak_detector: string;
    peak_laser: string;
    peak_color: string;
};

type NumericRow = {
    fluorophore: string;
    [key: string]: string | number;
};

type PanelPayload = {
    cytometer: string;
    libraries: LibraryInfo[];
    detectors: DetectorInfo[];
    fluorophores: FluorInfo[];
    selected: string[];
    spectra: NumericRow[];
    similarity: NumericRow[];
    complexity_index: number | null;
    peak_detectors: string[];
    error?: string;
};

type PanelExportResponse = {
    filename: string;
    content_type: string;
    content_base64: string;
    error?: string;
};

type TabId = 'panel' | 'similarity' | 'signatures';

type ImportedPanelRow = {
    fluor: string;
    marker: string;
};

const laserOrder = ['DeepUV', 'UV', 'Violet', 'Blue', 'YellowGreen', 'Red', 'IR', 'Other'];
const emptySlots = 18;

const formatMetric = (value: number | string | null | undefined) => {
    const numeric = typeof value === 'number' ? value : Number(value);
    if (!Number.isFinite(numeric)) return 'NA';
    return numeric.toFixed(2);
};

const toNumber = (value: string | number | undefined) => {
    if (typeof value === 'number') return Number.isFinite(value) ? value : 0;
    const parsed = Number(value);
    return Number.isFinite(parsed) ? parsed : 0;
};

const laserLabel = (laser: string) => {
    if (laser === 'YellowGreen') return 'YG 561';
    if (laser === 'Violet') return 'V 405';
    if (laser === 'Blue') return 'B 488';
    if (laser === 'Red') return 'R 640';
    if (laser === 'DeepUV') return 'DUV 320';
    if (laser === 'UV') return 'UV 355';
    if (laser === 'IR') return 'IR 781';
    return laser;
};

const unique = <T,>(values: T[]) => Array.from(new Set(values));

const linePath = (row: NumericRow, detectors: DetectorInfo[], width: number, height: number, left = 0) => {
    if (detectors.length === 0) return '';
    return detectors.map((det, index) => {
        const x = left + (detectors.length === 1 ? width / 2 : (index / (detectors.length - 1)) * width);
        const y = height - toNumber(row[det.detector]) * (height - 18) - 10;
        return `${index === 0 ? 'M' : 'L'}${x.toFixed(1)},${y.toFixed(1)}`;
    }).join(' ');
};

const similarityColor = (value: number) => {
    const alpha = Math.max(0.08, Math.min(0.82, value));
    const blue = Math.round(245 - alpha * 120);
    return `rgb(${blue}, ${Math.round(248 - alpha * 90)}, ${Math.round(252 - alpha * 35)})`;
};

const bandColor = (value: number) => {
    const v = Math.max(0, Math.min(1, value));
    const stops = [
        { at: 0, color: [0, 0, 255] },
        { at: 0.25, color: [0, 255, 255] },
        { at: 0.5, color: [0, 255, 0] },
        { at: 0.75, color: [255, 255, 0] },
        { at: 1, color: [255, 0, 0] },
    ];
    const upper = stops.findIndex(stop => v <= stop.at);
    const hi = stops[Math.max(1, upper)];
    const lo = stops[Math.max(0, Math.max(1, upper) - 1)];
    const range = Math.max(0.0001, hi.at - lo.at);
    const t = (v - lo.at) / range;
    const rgb = lo.color.map((channel, i) => Math.round(channel + (hi.color[i] - channel) * t));
    return `rgb(${rgb.join(', ')})`;
};

const signatureLogFromValue = (value: number) => {
    const v = Math.max(0, Math.min(1, value));
    return 0.35 + Math.pow(v, 0.72) * 5.65;
};

const signatureY = (logValue: number, top: number, height: number) => {
    const yPower = 1.5;
    const maxTransformed = Math.pow(6.35, yPower);
    const transformed = Math.pow(Math.max(0, Math.min(6.35, logValue)), yPower);
    return top + height - (transformed / maxTransformed) * height;
};

const signatureBandBins = (value: number) => {
    const offsets = [-0.36, -0.30, -0.24, -0.18, -0.12, -0.06, 0, 0.06, 0.12, 0.18, 0.24, 0.30, 0.36];
    const v = Math.max(0, Math.min(1, value));
    const center = signatureLogFromValue(v);
    return offsets.map((offset, index) => {
        const distance = Math.abs(index - (offsets.length - 1) / 2) / ((offsets.length - 1) / 2);
        const density = Math.max(0.05, 1 - distance) * Math.max(0.08, v);
        return {
            logValue: Math.max(0, Math.min(6, center + offset)),
            density,
        };
    });
};

const csvEscape = (value: string) => `"${value.replace(/"/g, '""')}"`;

const downloadBlob = (filename: string, blob: Blob) => {
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = filename;
    document.body.appendChild(a);
    a.click();
    a.remove();
    URL.revokeObjectURL(url);
};

const base64ToBlob = (content: string, contentType: string) => {
    const binary = window.atob(content);
    const bytes = new Uint8Array(binary.length);
    for (let i = 0; i < binary.length; i += 1) bytes[i] = binary.charCodeAt(i);
    return new Blob([bytes], { type: contentType });
};

const normalizeImportToken = (value: string) => value.toLowerCase().replace(/[^a-z0-9]+/g, '');

const parseDelimitedRows = (text: string, delimiter: string) => {
    const rows: string[][] = [];
    let row: string[] = [];
    let cell = '';
    let quoted = false;

    for (let i = 0; i < text.length; i += 1) {
        const char = text[i];
        const next = text[i + 1];

        if (char === '"') {
            if (quoted && next === '"') {
                cell += '"';
                i += 1;
            } else {
                quoted = !quoted;
            }
        } else if (char === delimiter && !quoted) {
            row.push(cell.trim());
            cell = '';
        } else if ((char === '\n' || char === '\r') && !quoted) {
            row.push(cell.trim());
            cell = '';
            if (row.some(value => value.length > 0)) rows.push(row);
            row = [];
            if (char === '\r' && next === '\n') i += 1;
        } else {
            cell += char;
        }
    }

    row.push(cell.trim());
    if (row.some(value => value.length > 0)) rows.push(row);
    if (rows.length > 0 && rows[0].length > 0) rows[0][0] = rows[0][0].replace(/^\uFEFF/, '');
    return rows;
};

const parseCsvLikeRows = (text: string) => {
    const delimiters = [',', '\t', ';'];
    const parsed = delimiters.map(delimiter => {
        const rows = parseDelimitedRows(text, delimiter);
        const multiColumnRows = rows.filter(row => row.length > 1).length;
        const totalCells = rows.reduce((sum, row) => sum + row.length, 0);
        return { delimiter, rows, score: multiColumnRows * 1000 + totalCells };
    });
    parsed.sort((a, b) => b.score - a.score);
    return parsed[0]?.rows || [];
};

const rowHasHeaderWords = (row: string[]) => row.some(value => {
    const key = normalizeImportToken(value);
    return ['marker', 'markers', 'antigen', 'target', 'targets', 'fluor', 'fluorophore', 'fluorochrome', 'dye', 'tag', 'color', 'colour', 'reagent'].includes(key);
});

const buildFluorLookup = (fluorophores: FluorInfo[]) => {
    const lookup = new Map<string, string>();
    fluorophores.forEach(info => {
        const canonical = info.fluorophore;
        const keys = [
            canonical,
            canonical.replace(/\s*\([^)]*\)/g, ''),
            canonical.replace(/^(.+?)\s*-\s*/g, '$1-'),
        ].map(normalizeImportToken).filter(Boolean);
        keys.forEach(key => {
            if (!lookup.has(key)) lookup.set(key, canonical);
        });
    });
    return lookup;
};

const matchImportedFluor = (value: string, lookup: Map<string, string>) => {
    const raw = value.trim();
    if (!raw) return '';
    const candidates = [
        raw,
        raw.replace(/\s*\([^)]*\)/g, ''),
        raw.split(/[;,|]/)[0] || raw,
    ].map(normalizeImportToken).filter(Boolean);

    for (const key of candidates) {
        const hit = lookup.get(key);
        if (hit) return hit;
    }
    return '';
};

const detectImportedPanelRows = (text: string, fluorophores: FluorInfo[]) => {
    const rows = parseCsvLikeRows(text);
    if (rows.length === 0) throw new Error('The imported CSV file is empty.');

    const lookup = buildFluorLookup(fluorophores);
    const firstRow = rows[0] || [];
    const firstRowHasFluor = firstRow.some(value => !!matchImportedFluor(value, lookup));
    const hasHeader = rowHasHeaderWords(firstRow) || (!firstRowHasFluor && rows.length > 1);
    const headers = hasHeader ? firstRow.map(normalizeImportToken) : [];
    const dataRows = hasHeader ? rows.slice(1) : rows;
    const maxCols = Math.max(...rows.map(row => row.length));

    const fluorScores = Array.from({ length: maxCols }, (_, colIndex) => {
        const valueMatches = dataRows.reduce((count, row) => count + (matchImportedFluor(row[colIndex] || '', lookup) ? 1 : 0), 0);
        const header = headers[colIndex] || '';
        const headerBonus = ['fluor', 'fluorophore', 'fluorochrome', 'dye', 'tag', 'color', 'colour', 'reagent'].includes(header) ? 2 : 0;
        return valueMatches + headerBonus;
    });
    const fluorCol = fluorScores.indexOf(Math.max(...fluorScores));
    if (fluorCol < 0 || fluorScores[fluorCol] <= 0) {
        throw new Error('No known fluorophores were recognized in the imported CSV for the selected cytometer.');
    }

    const markerHeaderHits = headers
        .map((header, index) => ({ header, index }))
        .filter(item => ['marker', 'markers', 'antigen', 'target', 'targets'].includes(item.header) && item.index !== fluorCol);
    const markerCol = markerHeaderHits[0]?.index ?? (maxCols > 1
        ? Array.from({ length: maxCols }, (_, index) => index).find(index => index !== fluorCol && fluorScores[index] === 0)
        : undefined);

    const imported: ImportedPanelRow[] = [];
    const seen = new Set<string>();
    dataRows.forEach(row => {
        const fluor = matchImportedFluor(row[fluorCol] || '', lookup);
        if (!fluor || seen.has(fluor)) return;
        seen.add(fluor);
        imported.push({
            fluor,
            marker: markerCol === undefined ? '' : (row[markerCol] || '').trim(),
        });
    });

    if (imported.length === 0) {
        throw new Error('No known fluorophores were recognized in the imported CSV for the selected cytometer.');
    }
    return imported;
};

const PanelBuilder = () => {
    const [payload, setPayload] = useState<PanelPayload | null>(null);
    const [cytometer, setCytometer] = useState('aurora');
    const [slots, setSlots] = useState<string[]>(Array(emptySlots).fill(''));
    const slotsRef = useRef<string[]>(Array(emptySlots).fill(''));
    const importInputRef = useRef<HTMLInputElement | null>(null);
    const [markers, setMarkers] = useState<Record<number, string>>({});
    const [queries, setQueries] = useState<Record<number, string>>({});
    const [activeSlot, setActiveSlot] = useState<number | null>(null);
    const [tab, setTab] = useState<TabId>('panel');
    const [loading, setLoading] = useState(true);
    const [error, setError] = useState('');
    const [exporting, setExporting] = useState(false);
    const [importing, setImporting] = useState(false);

    const selected = useMemo(() => slots.filter(Boolean), [slots]);

    useEffect(() => {
        slotsRef.current = slots;
    }, [slots]);

    const selectedSet = useMemo(() => new Set(selected), [selected]);

    const colorByFluor = useMemo(() => {
        const map = new Map<string, string>();
        payload?.fluorophores.forEach(f => map.set(f.fluorophore, f.peak_color || '#2688e8'));
        return map;
    }, [payload]);

    const fluorByName = useMemo(() => {
        const map = new Map<string, FluorInfo>();
        payload?.fluorophores.forEach(f => map.set(f.fluorophore, f));
        return map;
    }, [payload]);

    const spectraByName = useMemo(() => {
        const map = new Map<string, NumericRow>();
        payload?.spectra.forEach(row => map.set(row.fluorophore, row));
        return map;
    }, [payload]);

    const similarityByName = useMemo(() => {
        const map = new Map<string, NumericRow>();
        payload?.similarity.forEach(row => map.set(row.fluorophore, row));
        return map;
    }, [payload]);

    const lasers = useMemo(() => {
        if (!payload) return [];
        const present = unique(payload.detectors.map(d => d.laser));
        return laserOrder.filter(laser => present.includes(laser));
    }, [payload]);

    const emissions = useMemo(() => {
        if (!payload) return [];
        return unique(payload.detectors.map(d => d.emission)).sort((a, b) => a - b);
    }, [payload]);

    const selectedEntries = useMemo(() => {
        if (!payload) return [];
        return slots
            .map((fluor, slotIndex) => {
                const info = fluorByName.get(fluor);
                const peakDetector = info?.peak_detector || '';
                    const det = payload.detectors.find(d => d.detector === peakDetector);
                return {
                    fluor,
                    slotIndex,
                    marker: markers[slotIndex] || '',
                    color: info?.peak_color || '#2688e8',
                    peakLaser: det?.laser || info?.peak_laser || '',
                    peakEmission: det?.emission || 0,
                };
            })
            .filter(entry => entry.fluor);
    }, [fluorByName, markers, payload, slots]);

    const selectedRows = useMemo(() => slots
        .map((fluor, slotIndex) => ({
            fluor,
            marker: (markers[slotIndex] || '').trim(),
        }))
        .filter(row => row.fluor), [markers, slots]);

    const fetchPanel = async (nextCytometer: string, nextSelected: string[]) => {
        setError('');
        const res = await axios.post(`${API_BASE}/spectral_panel_metrics`, {
            cytometer: nextCytometer,
            fluorophores: nextSelected,
        });
        if (res.data?.error) throw new Error(String(res.data.error));
        const nextPayload = res.data as PanelPayload;
        setPayload(nextPayload);
        setCytometer(nextPayload.cytometer);
        return nextPayload;
    };

    useEffect(() => {
        const boot = async () => {
            try {
                const params = new URLSearchParams(window.location.search);
                const cytFromUrl = params.get('cytometer');
                const res = await axios.get(`${API_BASE}/spectral_panel${cytFromUrl ? `?cytometer=${encodeURIComponent(cytFromUrl)}` : ''}`);
                if (res.data?.error) throw new Error(String(res.data.error));
                const initial = res.data as PanelPayload;
                setPayload(initial);
                setCytometer(initial.cytometer);
                const initialSlots = [...initial.selected];
                while (initialSlots.length < emptySlots) initialSlots.push('');
                setSlots(initialSlots);
                slotsRef.current = initialSlots;
            } catch (err) {
                setError(err instanceof Error ? err.message : 'Could not load spectral libraries.');
            } finally {
                setLoading(false);
            }
        };
        void boot();
    }, []);

    const updateSlot = async (index: number, fluor: string) => {
        const currentSlots = slotsRef.current;
        if (fluor && currentSlots.some((existing, i) => i !== index && existing === fluor)) return;
        const nextSlots = currentSlots.map((existing, i) => (i === index ? fluor : existing));
        slotsRef.current = nextSlots;
        setSlots(nextSlots);
        setQueries(prev => ({ ...prev, [index]: '' }));
        setActiveSlot(null);
        await fetchPanel(cytometer, nextSlots.filter(Boolean)).catch(err => {
            setError(err instanceof Error ? err.message : 'Could not update panel.');
        });
    };

    const clearSlot = async (index: number) => {
        await updateSlot(index, '');
        setMarkers(prev => {
            const next = { ...prev };
            delete next[index];
            return next;
        });
    };

    const addSlot = () => setSlots(prev => [...prev, '']);

    const changeCytometer = async (nextCytometer: string) => {
        setLoading(true);
        try {
            const nextPayload = await fetchPanel(nextCytometer, []);
            const nextSlots = Array(emptySlots).fill('');
            slotsRef.current = nextSlots;
            setSlots(nextSlots);
            setMarkers({});
            setPayload(nextPayload);
        } catch (err) {
            setError(err instanceof Error ? err.message : 'Could not switch cytometer.');
        } finally {
            setLoading(false);
        }
    };

    const filteredOptions = (slotIndex: number) => {
        if (!payload) return [];
        const query = (queries[slotIndex] ?? '').trim().toLowerCase();
        return payload.fluorophores
            .filter(f => !selectedSet.has(f.fluorophore) || slots[slotIndex] === f.fluorophore)
            .filter(f => !query || f.fluorophore.toLowerCase().includes(query))
            .slice(0, 80);
    };

    const exportPanelCsv = () => {
        if (selectedRows.length === 0) {
            setError('Select at least one fluorophore before exporting a panel.');
            return;
        }
        const lines = [
            ['Marker', 'Fluorophore'].map(csvEscape).join(','),
            ...selectedRows.map(row => [row.marker, row.fluor].map(csvEscape).join(',')),
        ];
        downloadBlob(
            `spectreasy_${cytometer}_panel.csv`,
            new Blob([`${lines.join('\n')}\n`], { type: 'text/csv;charset=utf-8' })
        );
    };

    const exportPanelOverview = async () => {
        if (selectedRows.length === 0) {
            setError('Select at least one fluorophore before exporting a panel overview.');
            return;
        }
        setError('');
        setExporting(true);
        try {
            const res = await axios.post(`${API_BASE}/export_spectral_panel_overview`, {
                cytometer,
                fluorophores: selectedRows.map(row => row.fluor),
                markers: selectedRows.map(row => row.marker),
            });
            if (res.data?.error) throw new Error(String(res.data.error));
            const out = res.data as PanelExportResponse;
            downloadBlob(out.filename || `spectreasy_${cytometer}_panel_overview.pdf`, base64ToBlob(out.content_base64, out.content_type || 'application/pdf'));
        } catch (err) {
            setError(err instanceof Error ? err.message : 'Could not export panel overview.');
        } finally {
            setExporting(false);
        }
    };

    const importPanelCsv = async (file: File | null) => {
        if (!file || !payload) return;
        setError('');
        setImporting(true);
        try {
            const text = await file.text();
            const imported = detectImportedPanelRows(text, payload.fluorophores);
            const nextSlots = imported.map(row => row.fluor);
            while (nextSlots.length < emptySlots) nextSlots.push('');
            slotsRef.current = nextSlots;
            setSlots(nextSlots);
            const nextMarkers: Record<number, string> = {};
            imported.forEach((row, index) => {
                if (row.marker) nextMarkers[index] = row.marker;
            });
            setMarkers(nextMarkers);
            setQueries({});
            setActiveSlot(null);
            await fetchPanel(cytometer, imported.map(row => row.fluor));
        } catch (err) {
            setError(err instanceof Error ? err.message : 'Could not import panel CSV.');
        } finally {
            setImporting(false);
            if (importInputRef.current) importInputRef.current.value = '';
        }
    };

    if (loading) {
        return <div className="panel-builder"><div className="empty-state">Loading spectral panel builder...</div></div>;
    }

    if (error && !payload) {
        return <div className="panel-builder"><div className="error-state">{error}</div></div>;
    }

    if (!payload) return null;

    const chartWidth = Math.max(1040, payload.detectors.length * 22);
    const chartHeight = 230;
    const signatureLeft = 58;
    const signatureTop = 22;
    const signaturePlotHeight = 265;
    const signatureAxisBottom = signatureTop + signaturePlotHeight;
    const signatureHeight = signatureAxisBottom + 82;
    const signaturePlotWidth = chartWidth - signatureLeft - 18;
    const signatureColumnWidth = signaturePlotWidth / Math.max(1, payload.detectors.length);

    return (
        <div className="panel-builder">
            <header className="panel-topbar">
                <div>
                    <h1>Spectral Panel Builder</h1>
                    <p>{selected.length} fluorophores selected</p>
                </div>
                <div className="panel-actions">
                    <input
                        ref={importInputRef}
                        type="file"
                        accept=".csv,.tsv,.txt,text/csv,text/tab-separated-values"
                        className="hidden-file-input"
                        onChange={event => void importPanelCsv(event.target.files?.[0] || null)}
                    />
                    <button type="button" className="export-button" onClick={() => importInputRef.current?.click()} disabled={importing}>
                        <Upload size={16} />
                        {importing ? 'Importing...' : 'Import panel CSV'}
                    </button>
                    <button type="button" className="export-button" onClick={exportPanelCsv}>
                        <FileText size={16} />
                        Export panel CSV
                    </button>
                    <button type="button" className="export-button primary" onClick={() => void exportPanelOverview()} disabled={exporting}>
                        <Download size={16} />
                        {exporting ? 'Exporting...' : 'Export overview PDF'}
                    </button>
                </div>
            </header>
            <div className="panel-shell">
                <aside className="panel-sidebar">
                    <div className="panel-sidebar-head">
                        <select
                            className="instrument-select"
                            value={cytometer}
                            onChange={event => void changeCytometer(event.target.value)}
                        >
                            {payload.libraries.map(lib => (
                                <option key={lib.id} value={lib.id}>{lib.label}</option>
                            ))}
                        </select>
                    </div>
                    <div className="selector-list">
                        {slots.map((fluor, index) => {
                            const info = fluorByName.get(fluor);
                            const display = activeSlot === index ? (queries[index] ?? fluor) : fluor;
                            const color = info?.peak_color || '#d1d5db';
                            return (
                                <div className="selector-row" key={`slot-${index}`}>
                                    <div className="selector-swatch" style={{ background: color }} />
                                    <div>
                                        <input
                                            className="selector-input"
                                            value={display}
                                            placeholder="Select fluorophore"
                                            onFocus={() => {
                                                setActiveSlot(index);
                                                setQueries(prev => ({ ...prev, [index]: fluor }));
                                            }}
                                            onChange={event => {
                                                setActiveSlot(index);
                                                setQueries(prev => ({ ...prev, [index]: event.target.value }));
                                            }}
                                            onKeyDown={event => {
                                                if (event.key !== 'Enter') return;
                                                const first = filteredOptions(index)[0];
                                                if (first) void updateSlot(index, first.fluorophore);
                                            }}
                                        />
                                        {activeSlot === index && (
                                            <div className="fluor-dropdown">
                                                {filteredOptions(index).map(option => (
                                                    <button
                                                        type="button"
                                                        className="fluor-option"
                                                        key={option.fluorophore}
                                                        onMouseDown={event => event.preventDefault()}
                                                        onClick={() => void updateSlot(index, option.fluorophore)}
                                                    >
                                                        <span className="fluor-option-swatch" style={{ background: option.peak_color || '#d1d5db' }} />
                                                        {option.fluorophore}
                                                    </button>
                                                ))}
                                            </div>
                                        )}
                                    </div>
                                    <button className="clear-slot" type="button" onClick={() => void clearSlot(index)} aria-label="Clear fluorophore">
                                        {fluor ? 'x' : ''}
                                    </button>
                                </div>
                            );
                        })}
                        <button type="button" className="fluor-option" onClick={addSlot}>
                            <span><Plus size={16} /> Add fluorophore row</span>
                        </button>
                    </div>
                </aside>

                <main className="main-panel">
                    <div className="top-spectrum">
                        <svg className="spectrum-svg" width="100%" height={chartHeight + 56} viewBox={`0 0 ${chartWidth} ${chartHeight + 56}`} role="img" aria-label="Combined spectral signatures">
                            {[0, 25, 50, 75, 100].map(tick => {
                                const y = chartHeight - (tick / 100) * (chartHeight - 18) - 10;
                                return (
                                    <g key={tick}>
                                        <line x1={42} y1={y} x2={chartWidth - 8} y2={y} stroke="#d7dbe0" strokeWidth={1} />
                                        <text x={28} y={y + 4} fontSize={12} textAnchor="end">{tick}</text>
                                    </g>
                                );
                            })}
                            {payload.detectors.map((det, index) => {
                                const x = 42 + (index / Math.max(1, payload.detectors.length - 1)) * (chartWidth - 56);
                                return (
                                    <g key={det.detector}>
                                        <line x1={x} y1={6} x2={x} y2={chartHeight - 10} stroke="#d7dbe0" strokeWidth={1} />
                                        <text x={x} y={chartHeight + 10} fontSize={11} textAnchor="end" transform={`rotate(-55 ${x} ${chartHeight + 10})`}>
                                            {det.label}
                                        </text>
                                    </g>
                                );
                            })}
                            <line x1={42} y1={chartHeight - 10} x2={chartWidth - 8} y2={chartHeight - 10} stroke="#111" strokeWidth={4} />
                            {selected.map(fluor => {
                                const row = spectraByName.get(fluor);
                                if (!row) return null;
                                return (
                                    <path
                                        key={fluor}
                                        d={linePath(row, payload.detectors, chartWidth - 56, chartHeight, 42)}
                                        transform="translate(0 0)"
                                        fill="none"
                                        stroke={colorByFluor.get(fluor) || '#2688e8'}
                                        strokeWidth={2.2}
                                        strokeLinejoin="round"
                                        strokeLinecap="round"
                                        opacity={0.92}
                                    />
                                );
                            })}
                        </svg>
                    </div>

                    <div className="tabs-bar">
                        <button className={`tab-button ${tab === 'panel' ? 'active' : ''}`} onClick={() => setTab('panel')}>PANEL MATRIX</button>
                        <button className={`tab-button ${tab === 'similarity' ? 'active' : ''}`} onClick={() => setTab('similarity')}>SIMILARITY MATRIX</button>
                        <button className={`tab-button ${tab === 'signatures' ? 'active' : ''}`} onClick={() => setTab('signatures')}>SIGNATURES</button>
                        <select className="coexpression-select" defaultValue="">
                            <option value="">Show Co-Expression</option>
                        </select>
                        <div className="complexity-badge">Complexity Index: {formatMetric(payload.complexity_index)}</div>
                    </div>

                    {error && <div className="error-state">{error}</div>}

                    <section className="tab-content">
                        {tab === 'panel' && (
                            <div className="panel-matrix-wrap">
                                <table className="panel-matrix">
                                    <thead>
                                        <tr>
                                            <th className="emission-head" />
                                            {lasers.map(laser => (
                                                <th className="laser-head" key={laser} colSpan={2} style={{ background: payload.detectors.find(d => d.laser === laser)?.color || '#64748b' }}>
                                                    {laserLabel(laser)}
                                                </th>
                                            ))}
                                        </tr>
                                        <tr>
                                            <th className="emission-head">Emission</th>
                                            {lasers.flatMap(laser => [
                                                <th key={`${laser}-marker`}>Marker</th>,
                                                <th key={`${laser}-fluor`}>Fluor</th>,
                                            ])}
                                        </tr>
                                    </thead>
                                    <tbody>
                                        {emissions.map(emission => (
                                            <tr key={emission}>
                                                <td className="emission-cell">{emission}</td>
                                                {lasers.flatMap(laser => {
                                                    const entries = selectedEntries.filter(entry => entry.peakLaser === laser && entry.peakEmission === emission);
                                                    const occupied = entries.length > 0;
                                                    return [
                                                        <td key={`${emission}-${laser}-marker`} className={occupied ? '' : 'empty-region'}>
                                                            {entries.map(entry => (
                                                                <input
                                                                    key={entry.slotIndex}
                                                                    className="matrix-marker-input"
                                                                    value={entry.marker}
                                                                    placeholder="Marker"
                                                                    onChange={event => setMarkers(prev => ({ ...prev, [entry.slotIndex]: event.target.value }))}
                                                                />
                                                            ))}
                                                        </td>,
                                                        <td key={`${emission}-${laser}-fluor`} className={occupied ? '' : 'empty-region'}>
                                                            {entries.map(entry => (
                                                                <span className="matrix-fluor" key={entry.slotIndex} style={{ color: entry.color }}>{entry.fluor}</span>
                                                            ))}
                                                        </td>,
                                                    ];
                                                })}
                                            </tr>
                                        ))}
                                    </tbody>
                                </table>
                            </div>
                        )}

                        {tab === 'similarity' && (
                            <div>
                                <p style={{ margin: '0 0 18px', fontWeight: 600 }}>Double click the similarity panel to expand.</p>
                                <div className="similarity-wrap">
                                    {selected.length === 0 ? (
                                        <div className="empty-state">Select fluorophores to calculate the similarity matrix.</div>
                                    ) : (
                                        <table className="similarity-table">
                                            <tbody>
                                                {selected.map((rowName, rowIndex) => (
                                                    <tr key={rowName}>
                                                        <th className="row-label">{rowName}</th>
                                                        {selected.map((colName, colIndex) => {
                                                            if (colIndex > rowIndex) return null;
                                                            const value = rowName === colName ? 1 : toNumber(similarityByName.get(rowName)?.[colName]);
                                                            return (
                                                                <td key={colName} style={{ background: rowName === colName ? '#d1d1d1' : similarityColor(value), color: value > 0.5 ? '#fff' : '#111' }}>
                                                                    {value.toFixed(value === 1 ? 0 : 2)}
                                                                </td>
                                                            );
                                                        })}
                                                    </tr>
                                                ))}
                                                <tr>
                                                    <th />
                                                    {selected.map(name => (
                                                        <th className="col-label" key={name}><span className="rotated-label">{name}</span></th>
                                                    ))}
                                                </tr>
                                            </tbody>
                                        </table>
                                    )}
                                </div>
                                <p style={{ margin: '12px 0 0 14px', color: '#5b8fc7', fontWeight: 800 }}>Complexity Index: {formatMetric(payload.complexity_index)}</p>
                            </div>
                        )}

                        {tab === 'signatures' && (
                            <div className="signatures-wrap">
                                {selected.length === 0 ? (
                                    <div className="empty-state">Select fluorophores to view signatures.</div>
                                ) : selected.map(fluor => {
                                    const row = spectraByName.get(fluor);
                                    if (!row) return null;
                                    return (
                                        <div className="signature-card" key={fluor}>
                                            <h3>{fluor}</h3>
                                            <svg className="signature-band-svg" width="100%" height={signatureHeight} viewBox={`0 0 ${chartWidth} ${signatureHeight}`} role="img" aria-label={`${fluor} signature`}>
                                                <rect x={signatureLeft} y={signatureTop} width={signaturePlotWidth} height={signaturePlotHeight} fill="#ffffff" stroke="#e5e7eb" />
                                                {[0, 1, 2, 3, 4, 5, 6].map(tick => {
                                                    const y = signatureY(tick, signatureTop, signaturePlotHeight);
                                                    return (
                                                        <g key={tick}>
                                                            <line x1={signatureLeft} y1={y} x2={signatureLeft + signaturePlotWidth} y2={y} stroke="#e8ebef" strokeWidth={1} />
                                                            <text x={signatureLeft - 9} y={y + 4} fontSize={11} textAnchor="end">{`10^${tick}`}</text>
                                                        </g>
                                                    );
                                                })}
                                                <text x={14} y={signatureTop + signaturePlotHeight / 2} fontSize={13} fontWeight={700} textAnchor="middle" transform={`rotate(-90 14 ${signatureTop + signaturePlotHeight / 2})`}>Intensity</text>
                                                {payload.detectors.map((det, index) => {
                                                    const x = signatureLeft + index * signatureColumnWidth;
                                                    const centerX = x + signatureColumnWidth / 2;
                                                    const value = toNumber(row[det.detector]);
                                                    return (
                                                        <g key={det.detector}>
                                                            <line x1={centerX} y1={signatureTop} x2={centerX} y2={signatureAxisBottom} stroke="#f1f3f6" strokeWidth={1} />
                                                            {signatureBandBins(value).map((bin, bandIndex) => {
                                                                const y = signatureY(bin.logValue, signatureTop, signaturePlotHeight);
                                                                return (
                                                                    <rect
                                                                        key={`${det.detector}-${bandIndex}`}
                                                                        x={centerX - Math.max(3, signatureColumnWidth * 0.28)}
                                                                        y={y - 2.3}
                                                                        width={Math.max(4, signatureColumnWidth * 0.56)}
                                                                        height={4.6}
                                                                        fill={bandColor(bin.density)}
                                                                        opacity={0.95}
                                                                    />
                                                                );
                                                            })}
                                                            <text x={centerX} y={signatureAxisBottom + 15} fontSize={10} textAnchor="end" transform={`rotate(-65 ${centerX} ${signatureAxisBottom + 15})`}>{det.label}</text>
                                                        </g>
                                                    );
                                                })}
                                            </svg>
                                        </div>
                                    );
                                })}
                            </div>
                        )}
                    </section>
                </main>
            </div>
        </div>
    );
};

export default PanelBuilder;
