import { useEffect, useMemo, useRef, useState } from 'react';
import type { CSSProperties, PointerEvent as ReactPointerEvent } from 'react';
import axios from 'axios';
import { Moon, PanelLeftClose, PanelLeftOpen, Plus, Save, Sun, Trash2, Upload } from 'lucide-react';
import './PanelBuilder.css';
import { resolveApiBase } from './apiBase';

const API_BASE = resolveApiBase();

const unboxGuiState = (value: unknown): unknown => {
    if (Array.isArray(value)) {
        if (value.length === 1) return unboxGuiState(value[0]);
        return value.map(unboxGuiState);
    }
    if (value && typeof value === 'object') {
        return Object.fromEntries(Object.entries(value).map(([key, item]) => [key, unboxGuiState(item)]));
    }
    return value;
};

const normalizeMarkers = (value: unknown): Record<number, string> => {
    if (!value || typeof value !== 'object' || Array.isArray(value)) return {};
    return Object.fromEntries(
        Object.entries(value).map(([key, marker]) => [Number(key), String(unboxGuiState(marker) ?? '')])
    );
};

const PdfIcon = ({ size = 20 }: { size?: number }) => (
    <svg width={size + 4} height={size} viewBox="0 0 30 24" fill="none" aria-hidden="true">
        <path d="M4 2.75h13l5 5v13.5H4z" stroke="currentColor" strokeWidth="1.5" strokeLinejoin="round" />
        <path d="M17 2.75v5h5" stroke="currentColor" strokeWidth="1.5" strokeLinejoin="round" />
        <text x="13" y="16.4" fill="currentColor" fontSize="6.2" fontWeight="800" textAnchor="middle" fontFamily="Avenir Next, sans-serif">PDF</text>
    </svg>
);

type LibraryInfo = {
    id: string;
    label: string;
};

type ConfigurationInfo = {
    id: string;
    label: string;
    description: string;
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
    configuration: string;
    libraries: LibraryInfo[];
    configurations: ConfigurationInfo[];
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

const toSimilarityValue = (value: string | number | undefined) => {
    const numeric = toNumber(value);
    if (Math.abs(numeric) < 0.005) return 0;
    return Math.max(0, Math.min(1, numeric));
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

const detectorPointX = (index: number, count: number, left: number, width: number) => (
    count <= 1 ? left + width / 2 : left + (index / (count - 1)) * width
);

const detectorColumnCenterX = (index: number, count: number, left: number, width: number) => (
    left + ((index + 0.5) / Math.max(1, count)) * width
);

const linePath = (row: NumericRow, detectors: DetectorInfo[], width: number, height: number, left = 0) => {
    if (detectors.length === 0) return '';
    return detectors.map((det, index) => {
        const x = detectorPointX(index, detectors.length, left, width);
        const y = height - toNumber(row[det.detector]) * (height - 32) - 24;
        return `${index === 0 ? 'M' : 'L'}${x.toFixed(1)},${y.toFixed(1)}`;
    }).join(' ');
};

const wavelengthToColor = (wavelength: number) => {
    let r = 0, g = 0, b = 0;
    if (wavelength >= 350 && wavelength < 440) {
        r = -(wavelength - 440) / (440 - 350);
        b = 1.0;
    } else if (wavelength >= 440 && wavelength < 490) {
        g = (wavelength - 440) / (490 - 440);
        b = 1.0;
    } else if (wavelength >= 490 && wavelength < 510) {
        g = 1.0;
        b = -(wavelength - 510) / (510 - 490);
    } else if (wavelength >= 510 && wavelength < 580) {
        r = (wavelength - 510) / (580 - 510);
        g = 1.0;
    } else if (wavelength >= 580 && wavelength < 645) {
        r = 1.0;
        g = -(wavelength - 645) / (645 - 580);
    } else if (wavelength >= 645 && wavelength <= 780) {
        r = 1.0;
    } else if (wavelength > 780) {
        r = 0.5;
        b = 0.2;
    }
    
    let factor = 1.0;
    if (wavelength >= 350 && wavelength < 420) {
        factor = 0.3 + 0.7 * (wavelength - 350) / (420 - 350);
    } else if (wavelength > 700 && wavelength <= 780) {
        factor = 0.3 + 0.7 * (780 - wavelength) / (780 - 700);
    } else if (wavelength > 780) {
        factor = 0.3;
    }
    
    return `rgb(${Math.round(r * factor * 255)}, ${Math.round(g * factor * 255)}, ${Math.round(b * factor * 255)})`;
};

const laserWavelength = (laser: string) => {
    const key = laser.toLowerCase();
    if (key === 'deepuv') return 320;
    if (key === 'uv') return 355;
    if (key === 'violet') return 405;
    if (key === 'blue') return 488;
    if (key === 'yellowgreen') return 561;
    if (key === 'red') return 640;
    if (key === 'ir') return 781;
    return 0;
};

const getSimilarityStyle = (value: number, isDiag: boolean, currentTheme: 'light' | 'dark') => {
    if (isDiag) {
        return {
            background: 'var(--bg-cell-diag)',
            color: currentTheme === 'dark' ? '#fff' : '#111'
        };
    }
    
    if (currentTheme === 'dark') {
        const opacity = Math.max(0.12, Math.min(0.92, value));
        const r = Math.round(15 + value * 124); // 15 -> 139
        const g = Math.round(23 + value * 69);  // 23 -> 92
        const b = Math.round(42 + value * 204); // 42 -> 246
        return {
            background: `rgba(${r}, ${g}, ${b}, ${opacity})`,
            color: '#fff',
            textShadow: value > 0.5 ? '0 0 4px rgba(255, 255, 255, 0.4)' : 'none'
        };
    } else {
        const alpha = Math.max(0.08, Math.min(0.82, value));
        const blue = Math.round(245 - alpha * 120);
        const rgbStr = `rgb(${blue}, ${Math.round(248 - alpha * 90)}, ${Math.round(252 - alpha * 35)})`;
        return {
            background: rgbStr,
            color: value > 0.4 ? '#fff' : '#1e293b'
        };
    }
};

const bandColor = (value: number) => {
    const v = Math.max(0, Math.min(1, Math.pow(value, 0.55)));
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

const getCytometerName = (val: unknown): string => {
    if (!val) return '';
    if (Array.isArray(val)) return String(val[0] || '');
    return String(val);
};

const binEmission = (val: number, cyt: string): number => {
    const c = getCytometerName(cyt).toLowerCase();
    if (c === 'id7000') {
        return Math.round(val / 20) * 20;
    }
    if (c === 'discover' || c === 'xenith') {
        return Math.round(val / 15) * 15;
    }
    return val;
};

const mapDetectorToEmission = (detectorName: string): number => {
    const clean = detectorName.replace(/-A$/, '').toUpperCase();
    if (clean.startsWith('UV')) {
        const idx = parseInt(clean.slice(2), 10);
        const uvMap: Record<number, number> = {
            1: 370,
            2: 395,
            3: 420,
            4: 440,
            5: 450,
            6: 480,
            7: 480,
            8: 500,
            9: 520,
            10: 550,
            11: 570,
            12: 580,
            13: 600,
            14: 660,
            15: 750,
            16: 800
        };
        return uvMap[idx] || 370;
    }
    if (clean.startsWith('V')) {
        const idx = parseInt(clean.slice(1), 10);
        const vMap: Record<number, number> = {
            1: 420,
            2: 440,
            3: 450,
            4: 480,
            5: 480,
            6: 500,
            7: 550,
            8: 570,
            9: 580,
            10: 600,
            11: 660,
            12: 680,
            13: 690,
            14: 700,
            15: 730,
            16: 780
        };
        return vMap[idx] || 420;
    }
    if (clean.startsWith('B')) {
        const idx = parseInt(clean.slice(1), 10);
        const bMap: Record<number, number> = {
            1: 500,
            2: 520,
            3: 550,
            4: 550,
            5: 570,
            6: 580,
            7: 600,
            8: 600,
            9: 660,
            10: 680,
            11: 690,
            12: 700,
            13: 750,
            14: 780
        };
        return bMap[idx] || 500;
    }
    if (clean.startsWith('YG')) {
        const idx = parseInt(clean.slice(2), 10);
        const ygMap: Record<number, number> = {
            1: 570,
            2: 580,
            3: 600,
            4: 600,
            5: 660,
            6: 680,
            7: 700,
            8: 730,
            9: 750,
            10: 780
        };
        return ygMap[idx] || 570;
    }
    if (clean.startsWith('R')) {
        const idx = parseInt(clean.slice(1), 10);
        const rMap: Record<number, number> = {
            1: 660,
            2: 680,
            3: 700,
            4: 730,
            5: 730,
            6: 750,
            7: 780,
            8: 800
        };
        return rMap[idx] || 660;
    }
    return 0;
};

const PanelBuilder = () => {
    const [payload, setPayload] = useState<PanelPayload | null>(null);
    const [cytometer, setCytometer] = useState(() => getCytometerName(localStorage.getItem('spectreasy_cytometer') || 'aurora'));
    const [configuration, setConfiguration] = useState(() => getCytometerName(localStorage.getItem('spectreasy_configuration') || '5l_uv_v_b_yg_r'));
    const [slots, setSlots] = useState<string[]>(() => {
        try {
            const stored = localStorage.getItem('spectreasy_slots');
            if (stored) {
                const parsed = JSON.parse(stored);
                if (Array.isArray(parsed)) return parsed;
            }
        } catch {
            return Array(emptySlots).fill('');
        }
        return Array(emptySlots).fill('');
    });
    const slotsRef = useRef<string[]>(slots);
    const importInputRef = useRef<HTMLInputElement | null>(null);
    const [markers, setMarkers] = useState<Record<number, string>>(() => {
        try {
            const stored = localStorage.getItem('spectreasy_markers');
            return stored ? normalizeMarkers(JSON.parse(stored)) : {};
        } catch {
            return {};
        }
    });
    const bootDefaultsRef = useRef({ cytometer, configuration, slots, markers });
    const [queries, setQueries] = useState<Record<number, string>>({});
    const [activeSlot, setActiveSlot] = useState<number | null>(null);
    const [tab, setTab] = useState<TabId>('panel');
    const [loading, setLoading] = useState(true);
    const [error, setError] = useState('');
    const [exporting, setExporting] = useState(false);
    const [importing, setImporting] = useState(false);
    const [hoveredFluor, setHoveredFluor] = useState<string | null>(null);
    const [theme, setTheme] = useState<'light' | 'dark'>(() => {
        const stored = localStorage.getItem('spectreasy-theme') || localStorage.getItem('spectreasy_theme');
        if (stored === 'light' || stored === 'dark') return stored;
        return window.matchMedia('(prefers-color-scheme: dark)').matches ? 'dark' : 'light';
    });
    const [guiStateLoaded, setGuiStateLoaded] = useState(false);
    const [sidebarWidth, setSidebarWidth] = useState(214);
    const [sidebarCollapsed, setSidebarCollapsed] = useState(false);
    const [showPdfConfirm, setShowPdfConfirm] = useState(false);

    useEffect(() => {
        localStorage.setItem('spectreasy_cytometer', getCytometerName(cytometer));
    }, [cytometer]);

    useEffect(() => {
        localStorage.setItem('spectreasy_configuration', getCytometerName(configuration));
    }, [configuration]);

    useEffect(() => {
        localStorage.setItem('spectreasy-theme', theme);
        localStorage.removeItem('spectreasy_theme');
        document.documentElement.dataset.theme = theme;
    }, [theme]);

    useEffect(() => {
        localStorage.setItem('spectreasy_slots', JSON.stringify(slots));
    }, [slots]);

    useEffect(() => {
        localStorage.setItem('spectreasy_markers', JSON.stringify(markers));
    }, [markers]);

    useEffect(() => {
        if (!guiStateLoaded) return;
        const timer = window.setTimeout(() => {
            void axios.post(`${API_BASE}/gui_state`, {
                module: 'panel_builder',
                config_json: {
                    cytometer: getCytometerName(cytometer),
                    configuration: getCytometerName(configuration),
                    theme,
                    slots,
                    markers,
                    tab,
                    sidebarWidth,
                    sidebarCollapsed
                }
            }).catch(() => null);
        }, 500);
        return () => window.clearTimeout(timer);
    }, [cytometer, configuration, theme, slots, markers, tab, sidebarWidth, sidebarCollapsed, guiStateLoaded]);

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
        const cytName = getCytometerName(cytometer).toLowerCase();
        if (cytName === 'aurora') {
            return [370, 395, 420, 440, 450, 480, 500, 520, 550, 570, 580, 600, 660, 680, 690, 700, 730, 750, 780, 800];
        }
        if (!payload) return [];
        const rawEmissions = payload.detectors.map(d => binEmission(d.emission, cytometer));
        return unique(rawEmissions).sort((a, b) => a - b);
    }, [payload, cytometer]);

    const selectedEntries = useMemo(() => {
        if (!payload) return [];
        return slots
            .map((fluor, slotIndex) => {
                const info = fluorByName.get(fluor);
                const peakDetector = info?.peak_detector || '';
                const det = payload.detectors.find(d => d.detector === peakDetector);
                const peakLaser = det?.laser || info?.peak_laser || '';
                let peakEmission = det?.emission || 0;
                const cytName = getCytometerName(cytometer).toLowerCase();
                if (cytName === 'aurora' && peakDetector) {
                    peakEmission = mapDetectorToEmission(peakDetector);
                } else if (peakEmission > 0) {
                    peakEmission = binEmission(peakEmission, cytometer);
                }
                return {
                    fluor,
                    slotIndex,
                    marker: markers[slotIndex] || '',
                    color: info?.peak_color || '#2688e8',
                    peakLaser,
                    peakEmission,
                };
            })
            .filter(entry => entry.fluor);
    }, [fluorByName, markers, payload, slots, cytometer]);

    const selectedRows = useMemo(() => slots
        .map((fluor, slotIndex) => ({
            fluor,
            marker: (markers[slotIndex] || '').trim(),
        }))
        .filter(row => row.fluor), [markers, slots]);

    const selectedConfigurationLabel = useMemo(() => {
        const hit = payload?.configurations.find(config => config.id === configuration);
        return hit?.label || '';
    }, [configuration, payload]);

    const fetchPanel = async (nextCytometer: string, nextConfiguration: string, nextSelected: string[]) => {
        setError('');
        const res = await axios.post(`${API_BASE}/spectral_panel_metrics`, {
            cytometer: nextCytometer,
            configuration: nextConfiguration,
            fluorophores: nextSelected,
        });
        if (res.data?.error) throw new Error(String(res.data.error));
        const nextPayload = res.data as PanelPayload;
        setPayload(nextPayload);
        setCytometer(getCytometerName(nextPayload.cytometer));
        setConfiguration(getCytometerName(nextPayload.configuration));
        return nextPayload;
    };

    useEffect(() => {
        const boot = async () => {
            try {
                const stateRes = await axios.get(`${API_BASE}/gui_state?module=panel_builder`).catch(() => null);
                const saved = unboxGuiState(stateRes?.data?.config || {}) as Record<string, unknown>;
                const defaults = bootDefaultsRef.current;
                const savedCytometer = typeof saved.cytometer === 'string' ? getCytometerName(saved.cytometer) : defaults.cytometer;
                const savedConfiguration = typeof saved.configuration === 'string' ? getCytometerName(saved.configuration) : defaults.configuration;
                const savedSlots = Array.isArray(saved.slots) ? saved.slots.map(String) : defaults.slots;
                const savedMarkers = saved.markers && typeof saved.markers === 'object' ? normalizeMarkers(saved.markers) : defaults.markers;
                if (saved.theme === 'light' || saved.theme === 'dark') setTheme(saved.theme);
                if (saved.tab === 'panel' || saved.tab === 'similarity' || saved.tab === 'signatures') setTab(saved.tab);
                if (typeof saved.sidebarWidth === 'number' && Number.isFinite(saved.sidebarWidth)) {
                    setSidebarWidth(Math.min(440, Math.max(180, saved.sidebarWidth)));
                }
                if (typeof saved.sidebarCollapsed === 'boolean') setSidebarCollapsed(saved.sidebarCollapsed);
                setSlots(savedSlots);
                slotsRef.current = savedSlots;
                setMarkers(savedMarkers);
                const res = await axios.post(`${API_BASE}/spectral_panel_metrics`, {
                    cytometer: savedCytometer,
                    configuration: savedConfiguration,
                    fluorophores: savedSlots.filter(Boolean),
                });
                if (res.data?.error) throw new Error(String(res.data.error));
                const initial = res.data as PanelPayload;
                setPayload(initial);
                setCytometer(getCytometerName(initial.cytometer));
                setConfiguration(getCytometerName(initial.configuration));
                setGuiStateLoaded(true);
            } catch (err) {
                setError(err instanceof Error ? err.message : 'Could not load spectral libraries.');
            } finally {
                setLoading(false);
            }
        };
        void boot();
    }, []);

    useEffect(() => {
        const handleClickOutside = (event: MouseEvent) => {
            const target = event.target as HTMLElement;
            if (target && typeof target.closest === 'function') {
                if (!target.closest('.selector-row') && !target.closest('.clear-slot')) {
                    setActiveSlot(null);
                }
            }
        };
        document.addEventListener('mousedown', handleClickOutside);
        return () => document.removeEventListener('mousedown', handleClickOutside);
    }, []);

    const updateSlot = async (index: number, fluor: string) => {
        const currentSlots = slotsRef.current;
        if (fluor && currentSlots.some((existing, i) => i !== index && existing === fluor)) return;
        const nextSlots = currentSlots.map((existing, i) => (i === index ? fluor : existing));
        slotsRef.current = nextSlots;
        setSlots(nextSlots);
        localStorage.setItem('spectreasy_slots', JSON.stringify(nextSlots));
        setQueries(prev => ({ ...prev, [index]: '' }));
        setActiveSlot(null);
        await fetchPanel(cytometer, configuration, nextSlots.filter(Boolean)).catch(err => {
            setError(err instanceof Error ? err.message : 'Could not update panel.');
        });
    };

    const clearSlot = async (index: number) => {
        await updateSlot(index, '');
        setMarkers(prev => {
            const next = { ...prev };
            delete next[index];
            localStorage.setItem('spectreasy_markers', JSON.stringify(next));
            return next;
        });
    };

    const addSlot = () => setSlots(prev => {
        const next = [...prev, ''];
        localStorage.setItem('spectreasy_slots', JSON.stringify(next));
        return next;
    });

    const changeCytometer = async (nextCytometer: string) => {
        setLoading(true);
        try {
            const activeSelected = slots.filter(Boolean);
            const nextPayload = await fetchPanel(nextCytometer, configuration, activeSelected);
            const availableSet = new Set(nextPayload.fluorophores.map(f => f.fluorophore));
            const nextSlots = slots.map(fluor => (availableSet.has(fluor) ? fluor : ''));
            const nextMarkers: Record<number, string> = {};
            Object.entries(markers).forEach(([key, val]) => {
                const idx = parseInt(key, 10);
                if (nextSlots[idx]) {
                    nextMarkers[idx] = val;
                }
            });
            slotsRef.current = nextSlots;
            setSlots(nextSlots);
            setMarkers(nextMarkers);
            localStorage.setItem('spectreasy_slots', JSON.stringify(nextSlots));
            localStorage.setItem('spectreasy_markers', JSON.stringify(nextMarkers));
            setPayload(nextPayload);
        } catch (err) {
            setError(err instanceof Error ? err.message : 'Could not switch cytometer.');
        } finally {
            setLoading(false);
        }
    };

    const applyPayloadAvailability = (nextPayload: PanelPayload) => {
        const availableSet = new Set(nextPayload.fluorophores.map(f => f.fluorophore));
        const nextSlots = slotsRef.current.map(fluor => (availableSet.has(fluor) ? fluor : ''));
        const nextMarkers: Record<number, string> = {};
        Object.entries(markers).forEach(([key, val]) => {
            const idx = parseInt(key, 10);
            if (nextSlots[idx]) nextMarkers[idx] = val;
        });
        slotsRef.current = nextSlots;
        setSlots(nextSlots);
        setMarkers(nextMarkers);
        localStorage.setItem('spectreasy_slots', JSON.stringify(nextSlots));
        localStorage.setItem('spectreasy_markers', JSON.stringify(nextMarkers));
        setPayload(nextPayload);
    };

    const changeConfiguration = async (nextConfiguration: string) => {
        setLoading(true);
        try {
            const nextPayload = await fetchPanel(cytometer, nextConfiguration, slotsRef.current.filter(Boolean));
            applyPayloadAvailability(nextPayload);
        } catch (err) {
            setError(err instanceof Error ? err.message : 'Could not switch Aurora configuration.');
        } finally {
            setLoading(false);
        }
    };

    const clearSelection = async () => {
        setError('');
        const emptySlotsArray = Array(emptySlots).fill('');
        slotsRef.current = emptySlotsArray;
        setSlots(emptySlotsArray);
        setMarkers({});
        localStorage.setItem('spectreasy_slots', JSON.stringify(emptySlotsArray));
        localStorage.setItem('spectreasy_markers', JSON.stringify({}));
        setQueries({});
        setActiveSlot(null);
        setLoading(true);
        try {
            await fetchPanel(cytometer, configuration, []);
        } catch (err) {
            setError(err instanceof Error ? err.message : 'Could not clear selection.');
        } finally {
            setLoading(false);
        }
    };

    const filteredOptions = (slotIndex: number) => {
        if (!payload) return [];
        const query = (queries[slotIndex] ?? '').trim().toLowerCase();
        const baseOptions = payload.fluorophores
            .filter(f => !selectedSet.has(f.fluorophore) || slots[slotIndex] === f.fluorophore);

        if (!query) {
            return baseOptions.slice(0, 80);
        }

        const score = (name: string): number => {
            const lower = name.toLowerCase();
            if (lower === query) return 1000;
            if (lower.startsWith(query)) {
                return 800 - lower.length;
            }
            if (lower.includes(' ' + query) || lower.includes('-' + query)) {
                return 600 - lower.length;
            }
            if (lower.includes(query)) {
                return 100 - lower.length;
            }
            return 0;
        };

        return baseOptions
            .filter(f => score(f.fluorophore) > 0)
            .sort((a, b) => score(b.fluorophore) - score(a.fluorophore))
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
            `spectreasy_${cytometer}_${configuration}_panel.csv`,
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
                configuration,
                fluorophores: selectedRows.map(row => row.fluor),
                markers: selectedRows.map(row => row.marker),
            });
            if (res.data?.error) throw new Error(String(res.data.error));
            const out = res.data as PanelExportResponse;
            downloadBlob(out.filename || `spectreasy_${cytometer}_${configuration}_panel_overview.pdf`, base64ToBlob(out.content_base64, out.content_type || 'application/pdf'));
        } catch (err) {
            setError(err instanceof Error ? err.message : 'Could not export panel overview.');
        } finally {
            setExporting(false);
        }
    };

    const beginSidebarResize = (event: ReactPointerEvent<HTMLDivElement>) => {
        if (sidebarCollapsed) return;
        event.preventDefault();
        const startX = event.clientX;
        const startWidth = sidebarWidth;
        const previousCursor = document.body.style.cursor;
        const previousUserSelect = document.body.style.userSelect;
        document.body.style.cursor = 'col-resize';
        document.body.style.userSelect = 'none';

        const move = (moveEvent: PointerEvent) => {
            setSidebarWidth(Math.min(440, Math.max(180, startWidth + moveEvent.clientX - startX)));
        };
        const finish = () => {
            window.removeEventListener('pointermove', move);
            window.removeEventListener('pointerup', finish);
            window.removeEventListener('pointercancel', finish);
            document.body.style.cursor = previousCursor;
            document.body.style.userSelect = previousUserSelect;
        };
        window.addEventListener('pointermove', move);
        window.addEventListener('pointerup', finish);
        window.addEventListener('pointercancel', finish);
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
            localStorage.setItem('spectreasy_slots', JSON.stringify(nextSlots));
            localStorage.setItem('spectreasy_markers', JSON.stringify(nextMarkers));
            setQueries({});
            setActiveSlot(null);
            await fetchPanel(cytometer, configuration, imported.map(row => row.fluor));
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
    const spectrumLeft = 42;
    const spectrumRight = chartWidth - 8;
    const spectrumPlotWidth = spectrumRight - spectrumLeft;
    const signatureLeft = 58;
    const signatureTop = 22;
    const signaturePlotHeight = 265;
    const signatureAxisBottom = signatureTop + signaturePlotHeight;
    const signatureHeight = signatureAxisBottom + 82;
    const signaturePlotWidth = chartWidth - signatureLeft - 18;
    const signatureColumnWidth = signaturePlotWidth / Math.max(1, payload.detectors.length);

    return (
        <div className={`panel-builder ${theme}`}>
            <header className="panel-topbar">
                <div>
                    <h1>Spectral Panel Builder</h1>
                    <p>{selected.length} fluorophores selected{selectedConfigurationLabel ? ` / ${selectedConfigurationLabel}` : ''}</p>
                </div>
                <div className="panel-actions">
                    <button 
                        type="button" 
                        className="export-button" 
                        onClick={() => setTheme(prev => prev === 'light' ? 'dark' : 'light')}
                        aria-label="Toggle theme"
                        style={{ padding: '0 10px', width: '40px' }}
                    >
                        {theme === 'light' ? <Moon size={16} /> : <Sun size={16} />}
                    </button>
                    <button
                        type="button"
                        className="export-button icon-only"
                        onClick={clearSelection}
                        disabled={selected.length === 0}
                        aria-label="Clear selection"
                        title="Clear selection"
                        style={{ color: 'var(--accent-2)' }}
                    >
                        <Trash2 size={16} />
                    </button>
                    <input
                        ref={importInputRef}
                        type="file"
                        accept=".csv,.tsv,.txt,text/csv,text/tab-separated-values"
                        className="hidden-file-input"
                        onChange={event => void importPanelCsv(event.target.files?.[0] || null)}
                    />
                    <button type="button" className="export-button icon-only" onClick={() => importInputRef.current?.click()} disabled={importing} aria-label={importing ? 'Importing panel CSV' : 'Import panel CSV'} title={importing ? 'Importing…' : 'Import panel CSV'}>
                        <Upload size={16} />
                    </button>
                    <button type="button" className="export-button icon-only" onClick={exportPanelCsv} aria-label="Export panel CSV" title="Export panel CSV">
                        <Save size={16} />
                    </button>
                    <button type="button" className="export-button primary icon-only" onClick={() => {
                        if (selectedRows.length === 0) {
                            void exportPanelOverview();
                            return;
                        }
                        setShowPdfConfirm(true);
                    }} disabled={exporting} aria-label={exporting ? 'Exporting overview PDF' : 'Export overview PDF'} title={exporting ? 'Exporting…' : 'Export overview PDF'}>
                        <PdfIcon size={20} />
                    </button>
                </div>
            </header>
            <div className="panel-shell">
                <aside
                    className={`panel-sidebar ${sidebarCollapsed ? 'is-collapsed' : ''}`}
                    style={{ '--panel-sidebar-width': `${sidebarWidth}px` } as CSSProperties}
                >
                    <button
                        type="button"
                        className="panel-sidebar-toggle"
                        onClick={() => setSidebarCollapsed(previous => !previous)}
                        aria-label={sidebarCollapsed ? 'Show fluorophore sidebar' : 'Hide fluorophore sidebar'}
                        title={sidebarCollapsed ? 'Show fluorophore sidebar' : 'Hide fluorophore sidebar'}
                    >
                        {sidebarCollapsed ? <PanelLeftOpen size={14} /> : <PanelLeftClose size={14} />}
                    </button>
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
                        {payload.configurations.length > 1 && (
                            <select
                                className="configuration-select"
                                value={configuration}
                                onChange={event => void changeConfiguration(event.target.value)}
                                aria-label="Detector configuration"
                            >
                                {payload.configurations.map(config => (
                                    <option key={config.id} value={config.id}>
                                        {config.label}
                                    </option>
                                ))}
                            </select>
                        )}
                    </div>
                    <div className="selector-list">
                        {slots.map((fluor, index) => {
                            const info = fluorByName.get(fluor);
                            const display = activeSlot === index ? (queries[index] ?? fluor) : fluor;
                            const color = info?.peak_color || '#d1d5db';
                            const isHovered = hoveredFluor === fluor;
                            return (
                                <div 
                                    className="selector-row" 
                                    key={`slot-${index}`}
                                    onMouseEnter={() => fluor && setHoveredFluor(fluor)}
                                    onMouseLeave={() => fluor && setHoveredFluor(null)}
                                    style={fluor && isHovered ? {
                                        borderColor: 'rgba(59, 130, 246, 0.5)',
                                        background: 'rgba(30, 41, 59, 0.35)',
                                        boxShadow: '0 4px 12px rgba(59, 130, 246, 0.15)'
                                    } : {}}
                                >
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

                {!sidebarCollapsed && (
                    <div
                        className="panel-sidebar-resizer"
                        role="separator"
                        aria-label="Resize fluorophore sidebar"
                        aria-orientation="vertical"
                        onPointerDown={beginSidebarResize}
                    />
                )}

                <main className="main-panel">
                    <div className="top-spectrum">
                        <svg className="spectrum-svg" width={chartWidth} height={chartHeight + 56} viewBox={`0 0 ${chartWidth} ${chartHeight + 56}`} role="img" aria-label="Combined spectral signatures">
                            <defs>
                                <linearGradient id="rainbow-axis-gradient" x1="0%" y1="0%" x2="100%" y2="0%">
                                    {payload.detectors.map((det, index) => {
                                        const pct = (index / Math.max(1, payload.detectors.length - 1)) * 100;
                                        return (
                                            <stop
                                                key={`stop-${det.detector}`}
                                                offset={`${pct.toFixed(2)}%`}
                                                stopColor={wavelengthToColor(det.emission)}
                                            />
                                        );
                                    })}
                                </linearGradient>
                            </defs>
                            {[0, 25, 50, 75, 100].map(tick => {
                                const y = chartHeight - (tick / 100) * (chartHeight - 32) - 24;
                                return (
                                    <g key={tick}>
                                        <line x1={spectrumLeft} y1={y} x2={spectrumRight} y2={y} stroke="var(--chart-grid)" strokeWidth={1} />
                                        <text x={28} y={y + 4} fontSize={12} textAnchor="end" className="chart-axis-text">{tick}</text>
                                    </g>
                                );
                            })}
                            {payload.detectors.map((det, index) => {
                                const x = detectorPointX(index, payload.detectors.length, spectrumLeft, spectrumPlotWidth);
                                return (
                                    <g key={det.detector}>
                                        <line x1={x} y1={6} x2={x} y2={chartHeight - 24} stroke="var(--chart-grid)" strokeWidth={1} />
                                        <text x={x} y={chartHeight + 4} fontSize={11} textAnchor="end" transform={`rotate(-90 ${x} ${chartHeight + 4})`} className="chart-axis-text">
                                            {det.label}
                                        </text>
                                    </g>
                                );
                            })}
                            <rect x={spectrumLeft} y={chartHeight - 14} width={spectrumPlotWidth} height={6} fill="url(#rainbow-axis-gradient)" />
                            <line x1={spectrumLeft} y1={chartHeight - 24} x2={spectrumRight} y2={chartHeight - 24} stroke="var(--chart-axis)" strokeWidth={3} />
                            {selected.map(fluor => {
                                const row = spectraByName.get(fluor);
                                if (!row) return null;
                                const isHovered = hoveredFluor === fluor;
                                const hasHoverActive = hoveredFluor !== null;
                                const pathData = linePath(row, payload.detectors, spectrumPlotWidth, chartHeight, spectrumLeft);
                                return (
                                    <g key={fluor}>
                                        <path
                                            d={pathData}
                                            fill="none"
                                            stroke={colorByFluor.get(fluor) || '#2688e8'}
                                            strokeWidth={isHovered ? 4.2 : 2.4}
                                            strokeLinejoin="round"
                                            strokeLinecap="round"
                                            opacity={hasHoverActive ? (isHovered ? 1.0 : 0.15) : 0.92}
                                            style={{
                                                transition: 'all 0.15s ease-in-out',
                                                pointerEvents: 'none',
                                                filter: isHovered ? `drop-shadow(0 0 5px ${colorByFluor.get(fluor) || '#2688e8'})` : 'none'
                                            }}
                                        />
                                        <path
                                            d={pathData}
                                            fill="none"
                                            stroke="transparent"
                                            strokeWidth={22}
                                            strokeLinejoin="round"
                                            strokeLinecap="round"
                                            onMouseEnter={() => setHoveredFluor(fluor)}
                                            onMouseLeave={() => setHoveredFluor(null)}
                                            style={{
                                                cursor: 'pointer',
                                                pointerEvents: 'stroke'
                                            }}
                                        />
                                    </g>
                                );
                            })}
                        </svg>
                    </div>

                    <div className="tabs-bar">
                        <button className={`tab-button ${tab === 'panel' ? 'active' : ''}`} onClick={() => setTab('panel')}>PANEL MATRIX</button>
                        <button className={`tab-button ${tab === 'similarity' ? 'active' : ''}`} onClick={() => setTab('similarity')}>SIMILARITY MATRIX</button>
                        <button className={`tab-button ${tab === 'signatures' ? 'active' : ''}`} onClick={() => setTab('signatures')}>SIGNATURES</button>
                        <div style={{ flex: 1 }} />
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
                                        {emissions.flatMap(emission => {
                                            const laserEntriesMap = new Map<string, typeof selectedEntries>();
                                            lasers.forEach(laser => {
                                                const entries = selectedEntries.filter(entry => entry.peakLaser === laser && entry.peakEmission === emission);
                                                laserEntriesMap.set(laser, entries);
                                            });
                                            const maxEntries = Math.max(1, ...lasers.map(laser => laserEntriesMap.get(laser)?.length || 0));
                                            
                                            return Array.from({ length: maxEntries }).map((_, subIndex) => (
                                                <tr key={`${emission}-${subIndex}`}>
                                                    {subIndex === 0 ? (
                                                        <td className="emission-cell" rowSpan={maxEntries}>{emission}</td>
                                                    ) : null}
                                                    {lasers.flatMap(laser => {
                                                        const entries = laserEntriesMap.get(laser) || [];
                                                        const entry = entries[subIndex];
                                                        const occupied = !!entry;
                                                        const isImpossible = emission < laserWavelength(laser);
                                                        const cellClass = occupied ? '' : (isImpossible ? 'impossible-region' : 'empty-region');
                                                        
                                                        return [
                                                            <td key={`${emission}-${laser}-marker-${subIndex}`} className={cellClass}>
                                                                {!isImpossible && occupied && (
                                                                    <input
                                                                        key={entry.slotIndex}
                                                                        className="matrix-marker-input"
                                                                        value={entry.marker}
                                                                        placeholder="Marker"
                                                                        onChange={event => {
                                                                            const val = event.target.value;
                                                                            setMarkers(prev => {
                                                                                const next = { ...prev, [entry.slotIndex]: val };
                                                                                localStorage.setItem('spectreasy_markers', JSON.stringify(next));
                                                                                return next;
                                                                            });
                                                                        }}
                                                                    />
                                                                )}
                                                            </td>,
                                                            <td key={`${emission}-${laser}-fluor-${subIndex}`} className={cellClass}>
                                                                {!isImpossible && occupied && (
                                                                    <span className="matrix-fluor" key={entry.slotIndex} style={{ color: entry.color }}>{entry.fluor}</span>
                                                                )}
                                                            </td>,
                                                        ];
                                                    })}
                                                </tr>
                                            ));
                                        })}
                                    </tbody>
                                </table>
                            </div>
                        )}

                        {tab === 'similarity' && (
                            <div>
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
                                                            const value = rowName === colName ? 1 : toSimilarityValue(similarityByName.get(rowName)?.[colName]);
                                                            const cellStyle = getSimilarityStyle(value, rowName === colName, theme);
                                                            return (
                                                                <td key={colName} style={cellStyle}>
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
                            </div>
                        )}

                        {tab === 'signatures' && (
                            <div className="signatures-wrap">
                                {selected.length === 0 ? (
                                    <div className="empty-state">Select fluorophores to view signatures.</div>
                                ) : selected.map(fluor => {
                                    const row = spectraByName.get(fluor);
                                    if (!row) return null;
                                    
                                    const isDark = theme === 'dark';
                                    const plotBg = isDark ? '#0b1110' : '#f8f7f3';
                                    const plotStroke = isDark ? '#52615b' : '#c7c3ba';
                                    const textHeading = isDark ? '#f0f3f2' : '#17201d';
                                    const gridH = isDark ? 'rgba(169, 183, 177, 0.16)' : 'rgba(109, 117, 111, 0.14)';
                                    const gridV = isDark ? 'rgba(169, 183, 177, 0.1)' : 'rgba(109, 117, 111, 0.09)';
                                    
                                    return (
                                        <div className="signature-card" key={fluor}>
                                            <h3>{fluor}</h3>
                                            <svg className="signature-band-svg" width={chartWidth} height={signatureHeight} viewBox={`0 0 ${chartWidth} ${signatureHeight}`} role="img" aria-label={`${fluor} signature`}>
                                                <rect x={signatureLeft} y={signatureTop} width={signaturePlotWidth} height={signaturePlotHeight} fill={plotBg} stroke={plotStroke} />
                                                {[0, 1, 2, 3, 4, 5, 6].map(tick => {
                                                    const y = signatureY(tick, signatureTop, signaturePlotHeight);
                                                    return (
                                                        <g key={tick}>
                                                            <line x1={signatureLeft} y1={y} x2={signatureLeft + signaturePlotWidth} y2={y} stroke={gridH} strokeWidth={1} />
                                                            <text x={signatureLeft - 9} y={y + 4} fontSize={11} textAnchor="end" className="chart-axis-text">{`10^${tick}`}</text>
                                                        </g>
                                                    );
                                                })}
                                                <text x={14} y={signatureTop + signaturePlotHeight / 2} fontSize={13} fontWeight={700} textAnchor="middle" transform={`rotate(-90 14 ${signatureTop + signaturePlotHeight / 2})`} fill={textHeading}>Intensity</text>
                                                {payload.detectors.map((det, index) => {
                                                    const centerX = detectorColumnCenterX(index, payload.detectors.length, signatureLeft, signaturePlotWidth);
                                                    const value = toNumber(row[det.detector]);
                                                    return (
                                                        <g key={det.detector}>
                                                            <line x1={centerX} y1={signatureTop} x2={centerX} y2={signatureAxisBottom} stroke={gridV} strokeWidth={1} />
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
                                                            <text x={centerX} y={signatureAxisBottom + 12} fontSize={10} textAnchor="end" transform={`rotate(-90 ${centerX} ${signatureAxisBottom + 12})`} className="chart-axis-text">{det.label}</text>
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
            {showPdfConfirm && (
                <div className="panel-confirm-overlay" onMouseDown={() => setShowPdfConfirm(false)}>
                    <div className="panel-confirm-modal" role="dialog" aria-modal="true" aria-labelledby="panel-pdf-confirm-title" onMouseDown={event => event.stopPropagation()}>
                        <h2 id="panel-pdf-confirm-title">Generate PDF report for panel?</h2>
                        <p>The report will contain the current panel overview and selected fluorophores.</p>
                        <div className="panel-confirm-actions">
                            <button type="button" className="panel-confirm-cancel" onClick={() => setShowPdfConfirm(false)}>Cancel</button>
                            <button type="button" className="panel-confirm-submit" onClick={() => {
                                setShowPdfConfirm(false);
                                void exportPanelOverview();
                            }}>Generate PDF</button>
                        </div>
                    </div>
                </div>
            )}
        </div>
    );
};

export default PanelBuilder;
