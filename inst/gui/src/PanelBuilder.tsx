import { useEffect, useMemo, useState } from 'react';
import axios from 'axios';
import { Plus } from 'lucide-react';
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

type TabId = 'panel' | 'similarity' | 'signatures';

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

const PanelBuilder = () => {
    const [payload, setPayload] = useState<PanelPayload | null>(null);
    const [cytometer, setCytometer] = useState('aurora');
    const [slots, setSlots] = useState<string[]>(Array(emptySlots).fill(''));
    const [markers, setMarkers] = useState<Record<number, string>>({});
    const [queries, setQueries] = useState<Record<number, string>>({});
    const [activeSlot, setActiveSlot] = useState<number | null>(null);
    const [tab, setTab] = useState<TabId>('panel');
    const [loading, setLoading] = useState(true);
    const [error, setError] = useState('');

    const selected = useMemo(() => slots.filter(Boolean), [slots]);

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
            } catch (err) {
                setError(err instanceof Error ? err.message : 'Could not load spectral libraries.');
            } finally {
                setLoading(false);
            }
        };
        void boot();
    }, []);

    const updateSlot = async (index: number, fluor: string) => {
        if (fluor && slots.some((existing, i) => i !== index && existing === fluor)) return;
        const nextSlots = slots.map((existing, i) => (i === index ? fluor : existing));
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

    if (loading) {
        return <div className="panel-builder"><div className="empty-state">Loading spectral panel builder...</div></div>;
    }

    if (error && !payload) {
        return <div className="panel-builder"><div className="error-state">{error}</div></div>;
    }

    if (!payload) return null;

    const chartWidth = Math.max(1040, payload.detectors.length * 22);
    const chartHeight = 230;

    return (
        <div className="panel-builder">
            <div className="panel-shell">
                <aside className="panel-sidebar">
                    <div className="panel-sidebar-head">
                        <h1>Spectral Panel Builder</h1>
                        <p>{selected.length} fluorophores selected</p>
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
                        <svg width={chartWidth} height={chartHeight + 56} role="img" aria-label="Combined spectral signatures">
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
                        <button className="tab-button disabled">SIR</button>
                        <button className="tab-button disabled">SSM</button>
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
                                            <svg width={chartWidth} height={250} role="img" aria-label={`${fluor} signature`}>
                                                <line x1={42} y1={214} x2={chartWidth - 8} y2={214} stroke="#111" strokeWidth={2} />
                                                <line x1={42} y1={16} x2={42} y2={214} stroke="#111" strokeWidth={2} />
                                                {payload.detectors.map((det, index) => {
                                                    const x = 42 + (index / Math.max(1, payload.detectors.length - 1)) * (chartWidth - 56);
                                                    const value = toNumber(row[det.detector]);
                                                    const barHeight = Math.max(2, value * 170);
                                                    return (
                                                        <g key={det.detector}>
                                                            <rect x={x - 5} y={214 - barHeight} width={10} height={barHeight} fill={det.color} opacity={0.78} />
                                                            <text x={x} y={232} fontSize={10} textAnchor="end" transform={`rotate(-55 ${x} 232)`}>{det.label}</text>
                                                        </g>
                                                    );
                                                })}
                                                <path
                                                    d={linePath(row, payload.detectors, chartWidth - 56, 210, 42)}
                                                    fill="none"
                                                    stroke={colorByFluor.get(fluor) || '#111'}
                                                    strokeWidth={2}
                                                    strokeLinejoin="round"
                                                    strokeLinecap="round"
                                                />
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
