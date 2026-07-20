/* eslint-disable react-refresh/only-export-components -- shared panel primitives intentionally colocate the PDF glyph with pure helpers */
import { resolveApiBase } from './apiBase';

const API_BASE = resolveApiBase();
const panelProjectPath = () => window.sessionStorage.getItem('spectreasy-project-path') || '';
const panelProjectUrl = (path: string) => `${API_BASE}${path}${panelProjectPath() ? `${path.includes('?') ? '&' : '?'}project_path=${encodeURIComponent(panelProjectPath())}` : ''}`;
const panelProjectBody = <T extends Record<string, unknown>>(body: T) => ({ ...body, projectPath: panelProjectPath() });

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

export {
    API_BASE,
    PdfIcon,
    bandColor,
    base64ToBlob,
    binEmission,
    csvEscape,
    detectorColumnCenterX,
    detectorPointX,
    detectImportedPanelRows,
    downloadBlob,
    emptySlots,
    formatMetric,
    getCytometerName,
    getSimilarityStyle,
    laserLabel,
    laserOrder,
    laserWavelength,
    linePath,
    mapDetectorToEmission,
    normalizeMarkers,
    panelProjectBody,
    panelProjectUrl,
    signatureBandBins,
    signatureY,
    toNumber,
    toSimilarityValue,
    unique,
    unboxGuiState,
    wavelengthToColor,
};
export type {
    ConfigurationInfo,
    DetectorInfo,
    FluorInfo,
    ImportedPanelRow,
    LibraryInfo,
    NumericRow,
    PanelExportResponse,
    PanelPayload,
    TabId,
};
