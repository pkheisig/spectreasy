import { useState, useEffect, useRef } from 'react';
import { Settings2, Activity, Save, RefreshCw, FileText, Check, Sliders, Palette, Circle, Sun, Moon, Maximize2 } from 'lucide-react';
import axios from 'axios';
import ResidualPlot from './ResidualPlot';

const API_BASE = (() => {
    const envBase = (import.meta.env.VITE_API_BASE as string | undefined)?.trim();
    if (envBase) return envBase.replace(/\/$/, '');
    if (typeof window !== 'undefined') {
        if (window.location.port === '5174') return 'http://localhost:8000';
        return window.location.origin.replace(/\/$/, '');
    }
    return 'http://localhost:8000';
})();

interface MatrixRow {
    Marker: string;
    [key: string]: any;
}

const COLOR_PALETTES = {
    default: ['#60a5fa', '#f472b6', '#34d399', '#fbbf24', '#a78bfa', '#fb7185', '#38bdf8', '#4ade80', '#facc15', '#c084fc'],
    neon: ['#00ff87', '#ff00ff', '#00ffff', '#ff6b6b', '#ffd93d', '#6bcb77', '#4d96ff', '#ff9f45', '#ee6ff8', '#45b7d1'],
    pastel: ['#a8d8ea', '#aa96da', '#fcbad3', '#ffffd2', '#b5eaea', '#ffc8dd', '#bde0fe', '#cdb4db', '#fde4cf', '#caffbf'],
    warm: ['#ff6b35', '#f7c59f', '#efa00b', '#d65108', '#ff8243', '#ffb627', '#ff7b54', '#ffa07a', '#f4845f', '#e76f51'],
    cool: ['#5390d9', '#7400b8', '#6930c3', '#5e60ce', '#4ea8de', '#48bfe3', '#56cfe1', '#64dfdf', '#72efdd', '#80ffdb']
};

const SCATTER_COLORS = ['#3b82f6', '#ef4444', '#22c55e', '#f59e0b', '#8b5cf6', '#ec4899', '#14b8a6', '#f97316'];
const App = () => {
    const [matrices, setMatrices] = useState<string[]>([]);
    const [currentFile, setCurrentFile] = useState('refined_reference_matrix.csv');
    const [sampleFiles, setSampleFiles] = useState<string[]>([]);
    const [currentSample, setCurrentSample] = useState('');
    const [sampleImportStatus, setSampleImportStatus] = useState<'idle' | 'importing' | 'loaded' | 'error'>('idle');
    const [sampleImportMessage, setSampleImportMessage] = useState('');
    const [configFiles, setConfigFiles] = useState<string[]>([]);
    const [currentConfig, setCurrentConfig] = useState('gui_config.json');
    const [matrix, setMatrix] = useState<MatrixRow[]>([]);
    const [loading, setLoading] = useState(true);
    const [detectors, setDetectors] = useState<string[]>([]);
    const [detectorLabels, setDetectorLabels] = useState<string[]>([]);
    const [selectedMarkers, setSelectedMarkers] = useState<string[]>([]);
    const [unmixedData, setUnmixedData] = useState<any[]>([]);
    const [rawData, setRawData] = useState<any[]>([]);
    const [isDragging, setIsDragging] = useState(false);
    const [isUnmixingMatrix, setIsUnmixingMatrix] = useState(false);
    const [saveStatus, setSaveStatus] = useState<'idle' | 'saving' | 'saved'>('idle');
    const [configStatus, setConfigStatus] = useState<'idle' | 'saving' | 'saved' | 'loading'>('idle');

    const [lineWidth, setLineWidth] = useState(0.8);
    const [lineOpacity, setLineOpacity] = useState(0.85);
    const [colorPalette, setColorPalette] = useState<keyof typeof COLOR_PALETTES>('default');
    const [showControls, setShowControls] = useState(true);

    const [signatureHeight, setSignatureHeight] = useState(300);
    const [signatureDetWidth, setSignatureDetWidth] = useState(26);
    const [residualCellSize, setResidualCellSize] = useState(130);
    const [pointSize, setPointSize] = useState(1.5);
    const [pointOpacity, setPointOpacity] = useState(0.5);
    const [pointColor, setPointColor] = useState('#3b82f6');
    const [dragSensitivity, setDragSensitivity] = useState(0.05);

    const [theme, setTheme] = useState<'dark' | 'light'>('light');
    const [pageScroll, setPageScroll] = useState(true);
    const sampleFileInputRef = useRef<HTMLInputElement>(null);

    const fetchMatrices = async () => {
        const res = await axios.get(`${API_BASE}/matrices`);
        const list = Array.isArray(res.data) ? res.data : [];
        setMatrices(list);
        return list;
    };

    const fetchSamples = async () => {
        const res = await axios.get(`${API_BASE}/samples`);
        const list = Array.isArray(res.data) ? res.data : [];
        setSampleFiles(list);
        return list;
    };

    const fetchConfigs = async () => {
        const res = await axios.get(`${API_BASE}/configs`);
        const list = Array.isArray(res.data) ? res.data : [];
        setConfigFiles(list);
        if (list.length > 0 && !list.includes(currentConfig)) {
            setCurrentConfig(list[0]);
        }
        return list;
    };

    const fetchSampleData = async (sampleName: string, matrixData: any[], filenameForType: string, detNames: string[]) => {
        const q = sampleName && sampleName.length > 0
            ? `?sample_name=${encodeURIComponent(sampleName)}`
            : '';
        const resData = await axios.get(`${API_BASE}/data${q}`);
        if (!resData.data.error) {
            const raw = resData.data.raw_data;
            setRawData(raw);
            setDetectorLabels(resData.data.detector_labels || detNames);
            if (resData.data.sample_name) setCurrentSample(resData.data.sample_name);
            await runUnmix(matrixData, raw, filenameForType);
        }
    };

    const fetchData = async (filename = currentFile, sampleName = currentSample) => {
        setLoading(true);
        const resMatrix = await axios.get(`${API_BASE}/load_matrix?filename=${filename}`);
        if (!resMatrix.data.error) {
            const matrixData = resMatrix.data;
            setMatrix(matrixData);
            const isUnmix = filename.toLowerCase().includes('unmixing') || filename.toLowerCase().includes('w_');
            setIsUnmixingMatrix(isUnmix);
            const detNames = Object.keys(matrixData[0]).filter(k => k !== 'Marker');
            setDetectors(detNames);
            const allMarkers = matrixData.map((r: any) => r.Marker);
            setSelectedMarkers(allMarkers);
            const activeSample = sampleName && sampleName.length > 0 ? sampleName : (sampleFiles[0] || '');
            await fetchSampleData(activeSample, matrixData, filename, detNames);
        }
        setCurrentFile(filename);
        setLoading(false);
    };

    useEffect(() => {
        const init = async () => {
            const [mats, samples] = await Promise.all([fetchMatrices(), fetchSamples()]);
            await fetchConfigs();
            const firstMatrix = mats.includes(currentFile) ? currentFile : (mats[0] || currentFile);
            const firstSample = samples[0] || '';
            if (firstSample) setCurrentSample(firstSample);
            await fetchData(firstMatrix, firstSample);
        };
        init();
    }, []);

    const runUnmix = async (currentM: any[], currentRaw: any[], filename = currentFile) => {
        const M_obj: any = {};
        currentM.forEach(row => {
            M_obj[row.Marker] = { ...row };
            delete M_obj[row.Marker].Marker;
        });
        let useType = 'reference';
        if (typeof filename === 'string') {
            if (filename.toLowerCase().includes('unmixing') || filename.toLowerCase().includes('w_')) useType = 'unmixing';
        } else if (typeof filename === 'boolean') {
            useType = filename ? 'unmixing' : 'reference';
        }
        const res = await axios.post(`${API_BASE}/unmix`, {
            matrix_json: M_obj,
            raw_data_json: currentRaw,
            type: useType
        });
        setUnmixedData(res.data);
    };

    const handleResidualAdjust = (xMarker: string, yMarker: string, alpha: number) => {
        const newMatrix = matrix.map(row => {
            if (row.Marker === yMarker) {
                const rowX = matrix.find(r => r.Marker === xMarker);
                if (!rowX) return row;
                const newRow = { ...row };
                Object.keys(row).forEach(key => {
                    if (key !== 'Marker' && typeof row[key] === 'number') {
                        newRow[key] = Number(row[key]) + alpha * Number(rowX[key]);
                    }
                });
                return newRow;
            }
            return row;
        });
        setMatrix(newMatrix);
        runUnmix(newMatrix, rawData, isUnmixingMatrix as any);
    };

    const saveMatrix = async () => {
        setSaveStatus('saving');
        const newName = currentFile.replace('.csv', '_adjusted.csv');
        const result = await axios.post(`${API_BASE}/save_matrix`, {
            filename: newName,
            matrix_json: matrix
        }).catch(() => null);
        if (result) {
            setSaveStatus('saved');
            fetchMatrices();
            setTimeout(() => setSaveStatus('idle'), 2000);
        } else {
            setSaveStatus('idle');
        }
    };

    const applyConfig = (cfg: any) => {
        if (!cfg || typeof cfg !== 'object') return;
        if (typeof cfg.lineWidth === 'number') setLineWidth(cfg.lineWidth);
        if (typeof cfg.lineOpacity === 'number') setLineOpacity(cfg.lineOpacity);
        if (typeof cfg.colorPalette === 'string' && cfg.colorPalette in COLOR_PALETTES) {
            setColorPalette(cfg.colorPalette as keyof typeof COLOR_PALETTES);
        }
        if (typeof cfg.showControls === 'boolean') setShowControls(cfg.showControls);
        if (typeof cfg.signatureHeight === 'number') setSignatureHeight(cfg.signatureHeight);
        if (typeof cfg.signatureDetWidth === 'number') setSignatureDetWidth(cfg.signatureDetWidth);
        if (typeof cfg.residualCellSize === 'number') setResidualCellSize(cfg.residualCellSize);
        if (typeof cfg.pointSize === 'number') setPointSize(cfg.pointSize);
        if (typeof cfg.pointOpacity === 'number') setPointOpacity(cfg.pointOpacity);
        if (typeof cfg.pointColor === 'string') setPointColor(cfg.pointColor);
        if (typeof cfg.dragSensitivity === 'number') {
            setDragSensitivity(Math.max(0.01, Math.min(0.1, cfg.dragSensitivity)));
        }
    };

    const buildConfig = () => ({
        lineWidth,
        lineOpacity,
        colorPalette,
        showControls,
        signatureHeight,
        signatureDetWidth,
        residualCellSize,
        pointSize,
        pointOpacity,
        pointColor,
        dragSensitivity
    });

    const saveViewConfig = async () => {
        const name = currentConfig.trim().length > 0 ? currentConfig.trim() : 'gui_config.json';
        setConfigStatus('saving');
        const result = await axios.post(`${API_BASE}/save_config`, {
            filename: name,
            config_json: buildConfig()
        }).catch(() => null);
        if (result?.data?.success) {
            const savedName = result.data.filename || name;
            setCurrentConfig(savedName);
            await fetchConfigs();
            setConfigStatus('saved');
            setTimeout(() => setConfigStatus('idle'), 1500);
        } else {
            setConfigStatus('idle');
        }
    };

    const loadViewConfig = async (name = currentConfig) => {
        if (!name || name.trim().length === 0) return;
        setConfigStatus('loading');
        const result = await axios.get(`${API_BASE}/load_config?filename=${encodeURIComponent(name)}`).catch(() => null);
        if (result?.data && !result.data.error) {
            applyConfig(result.data);
            setCurrentConfig(name);
        }
        setConfigStatus('idle');
    };

    const handleSampleChange = async (sampleName: string) => {
        setCurrentSample(sampleName);
        if (!sampleName || matrix.length === 0) return;
        setLoading(true);
        await fetchSampleData(sampleName, matrix, currentFile, detectors);
        setLoading(false);
    };

    const readFileAsBase64 = (file: File): Promise<string> =>
        new Promise((resolve, reject) => {
            const reader = new FileReader();
            reader.onload = () => {
                const result = reader.result;
                if (typeof result !== 'string') {
                    reject(new Error('Could not read selected file.'));
                    return;
                }
                const payload = result.split(',')[1] || '';
                if (!payload) {
                    reject(new Error('Selected file is empty.'));
                    return;
                }
                resolve(payload);
            };
            reader.onerror = () => reject(new Error('Failed reading selected file.'));
            reader.readAsDataURL(file);
        });

    const triggerSampleFilePicker = () => {
        sampleFileInputRef.current?.click();
    };

    const handleSampleFilePicked = async (event: React.ChangeEvent<HTMLInputElement>) => {
        const file = event.target.files?.[0];
        event.target.value = '';
        if (!file) return;

        setSampleImportStatus('importing');
        setSampleImportMessage(`Importing ${file.name}...`);

        try {
            const contentBase64 = await readFileAsBase64(file);
            const response = await axios.post(`${API_BASE}/import_sample_content`, {
                filename: file.name,
                content_base64: contentBase64
            });

            if (response.data?.error) {
                throw new Error(String(response.data.error));
            }

            const importedName = String(response.data?.filename || file.name);
            await fetchSamples();
            setSampleImportStatus('loaded');
            setSampleImportMessage(`Loaded ${importedName}`);

            if (matrix.length > 0) {
                await handleSampleChange(importedName);
            } else {
                setCurrentSample(importedName);
            }
        } catch (err) {
            const msg = err instanceof Error ? err.message : 'Failed to import sample file.';
            setSampleImportStatus('error');
            setSampleImportMessage(msg);
        }
    };

    const svgRef = useRef<SVGSVGElement>(null);
    const handleMouseMove = (e: React.MouseEvent) => {
        if (!isDragging || selectedMarkers.length !== 1) return;
        const marker = selectedMarkers[0];
        const svg = svgRef.current;
        if (!svg) return;
        const rect = svg.getBoundingClientRect();
        const x = e.clientX - rect.left;
        const y = e.clientY - rect.top;
        const detIdx = Math.round((x / rect.width) * (detectors.length - 1));
        const detName = detectors[detIdx];
        const yNorm = y / rect.height;
        const newVal = Math.max(0, Math.min(1, 1 - (yNorm - 0.05) / 0.9));
        const nextMatrix = matrix.map(row => {
            if (row.Marker === marker) return { ...row, [detName]: newVal };
            return row;
        });
        setMatrix(nextMatrix);
        runUnmix(nextMatrix, rawData, isUnmixingMatrix as any);
    };

    const colors = COLOR_PALETTES[colorPalette];
    const markerNames = matrix.map(m => m.Marker);
    const chartWidth = detectors.length * signatureDetWidth;

    // iOS 26 Glassy Theme
    const glassyTheme = theme === 'dark' ? {
        // Gentle gradient background - softer than before
        bgGradient: 'linear-gradient(135deg, #1a1f35 0%, #0f172a 25%, #1e1b4b 50%, #1f2937 75%, #111827 100%)',
        // Glass panel styles
        glassBg: 'rgba(255, 255, 255, 0.03)',
        glassBorder: 'rgba(255, 255, 255, 0.08)',
        glassHighlight: 'rgba(255, 255, 255, 0.05)',
        glassBlur: 'blur(20px)',
        // Text colors
        text: '#f1f5f9',
        textMuted: 'rgba(148, 163, 184, 0.8)',
        textDim: 'rgba(203, 213, 225, 0.7)',
        // Accent colors - softer teal/cyan
        accent: '#38bdf8',
        accentGlow: 'rgba(56, 189, 248, 0.15)',
        // Input styles
        inputBg: 'rgba(15, 23, 42, 0.6)',
        // Grid
        gridLine: 'rgba(100, 116, 139, 0.15)',
    } : {
        bgGradient: 'linear-gradient(135deg, #f8fafc 0%, #e2e8f0 25%, #f1f5f9 50%, #e5e7eb 75%, #f9fafb 100%)',
        glassBg: 'rgba(255, 255, 255, 0.7)',
        glassBorder: 'rgba(0, 0, 0, 0.08)',
        glassHighlight: 'rgba(255, 255, 255, 0.9)',
        glassBlur: 'blur(20px)',
        text: '#0f172a',
        textMuted: 'rgba(71, 85, 105, 0.8)',
        textDim: 'rgba(100, 116, 139, 0.9)',
        accent: '#0ea5e9',
        accentGlow: 'rgba(14, 165, 233, 0.15)',
        inputBg: 'rgba(241, 245, 249, 0.9)',
        gridLine: 'rgba(148, 163, 184, 0.2)',
    };

    const g = glassyTheme;

    // Glassmorphism card style
    const glassCard = {
        background: g.glassBg,
        backdropFilter: g.glassBlur,
        WebkitBackdropFilter: g.glassBlur,
        border: `1px solid ${g.glassBorder}`,
        borderRadius: '16px',
        boxShadow: `0 8px 32px rgba(0, 0, 0, 0.12), inset 0 1px 0 ${g.glassHighlight}`,
    };

    const glassButton = {
        background: g.glassBg,
        backdropFilter: g.glassBlur,
        WebkitBackdropFilter: g.glassBlur,
        border: `1px solid ${g.glassBorder}`,
        borderRadius: '12px',
    };

    if (loading) return (
        <div style={{
            height: '100vh',
            width: '100vw',
            display: 'flex',
            alignItems: 'center',
            justifyContent: 'center',
            background: g.bgGradient,
            color: g.accent,
            fontFamily: '-apple-system, BlinkMacSystemFont, "SF Pro Display", "Segoe UI", sans-serif',
            fontWeight: 500
        }}>
            <RefreshCw className="animate-spin" style={{ marginRight: 12 }} size={24} /> Initializing Workspace...
        </div>
    );

    const lowerTriangleCells = selectedMarkers.length * (selectedMarkers.length - 1) / 2;

    return (
        <div style={{
            minHeight: '100vh',
            height: pageScroll ? 'auto' : '100vh',
            width: '100vw',
            display: 'flex',
            flexDirection: 'column',
            background: g.bgGradient,
            color: g.text,
            fontFamily: '-apple-system, BlinkMacSystemFont, "SF Pro Display", "Segoe UI", sans-serif',
            overflow: pageScroll ? 'auto' : 'hidden'
        }}>
            {/* Navbar */}
            <header style={{
                height: 64,
                flexShrink: 0,
                borderBottom: `1px solid ${g.glassBorder}`,
                display: 'flex',
                alignItems: 'center',
                padding: '0 20px',
                justifyContent: 'space-between',
                background: g.glassBg,
                backdropFilter: g.glassBlur,
                WebkitBackdropFilter: g.glassBlur,
                position: 'sticky',
                top: 0,
                zIndex: 50
            }}>
                <div style={{ display: 'flex', alignItems: 'center', gap: 14 }}>
                    <div style={{
                        width: 44,
                        height: 44,
                        background: `linear-gradient(135deg, ${g.accent}, #8b5cf6)`,
                        borderRadius: 12,
                        display: 'flex',
                        alignItems: 'center',
                        justifyContent: 'center',
                        boxShadow: `0 4px 16px ${g.accentGlow}`
                    }}>
                        <Activity color="white" size={22} />
                    </div>
                    <div>
                        <h1 style={{ fontWeight: 700, fontSize: 18, margin: 0, letterSpacing: '-0.02em' }}>spectreasy</h1>
                        <p style={{ fontSize: 10, color: g.textMuted, margin: 0, letterSpacing: '0.1em', textTransform: 'uppercase' }}>Interactive Spectral Tuner</p>
                    </div>
                </div>
                <div style={{ display: 'flex', alignItems: 'center', gap: 12 }}>
                    <button onClick={() => setPageScroll(!pageScroll)} style={{ ...glassButton, padding: 10, color: pageScroll ? g.accent : g.textDim, cursor: 'pointer' }}>
                        <Maximize2 size={18} />
                    </button>
                    <button onClick={() => setTheme(theme === 'dark' ? 'light' : 'dark')} style={{ ...glassButton, padding: 10, color: g.textDim, cursor: 'pointer' }}>
                        {theme === 'dark' ? <Sun size={18} /> : <Moon size={18} />}
                    </button>
                    <button onClick={() => setShowControls(!showControls)} style={{ ...glassButton, padding: 10, color: showControls ? g.accent : g.textDim, cursor: 'pointer' }}>
                        <Sliders size={18} />
                    </button>
                    <button onClick={() => fetchData()} style={{ ...glassButton, padding: 10, color: g.textDim, cursor: 'pointer' }}>
                        <RefreshCw size={18} />
                    </button>
                    <button onClick={saveMatrix} style={{
                        background: `linear-gradient(135deg, ${g.accent}, #8b5cf6)`,
                        border: 'none',
                        borderRadius: 12,
                        padding: '10px 20px',
                        color: 'white',
                        fontWeight: 600,
                        fontSize: 14,
                        cursor: 'pointer',
                        display: 'flex',
                        alignItems: 'center',
                        gap: 8,
                        boxShadow: `0 4px 16px ${g.accentGlow}`
                    }}>
                        <Save size={16} /> {saveStatus === 'saving' ? 'Saving...' : saveStatus === 'saved' ? 'Saved!' : 'Save'}
                    </button>
                </div>
            </header>

            <div style={{ flex: 1, display: 'flex', overflow: pageScroll ? 'visible' : 'hidden' }}>
                {/* Left Sidebar */}
                <aside style={{
                    width: 320,
                    flexShrink: 0,
                    borderRight: `1px solid ${g.glassBorder}`,
                    background: g.glassBg,
                    backdropFilter: g.glassBlur,
                    WebkitBackdropFilter: g.glassBlur,
                    display: 'flex',
                    flexDirection: 'column',
                    overflow: 'auto',
                    position: 'sticky',
                    top: 64,
                    height: 'calc(100vh - 64px)'
                }}>
                    <div style={{ padding: 16 }}>
                        <h3 style={{ fontSize: 10, fontWeight: 700, color: g.textMuted, letterSpacing: '0.15em', textTransform: 'uppercase', marginBottom: 10 }}>Matrices</h3>
                        <div style={{ display: 'flex', flexDirection: 'column', gap: 4 }}>
                            {matrices.map(m => (
                                <button key={m} onClick={() => fetchData(m, currentSample)} style={{
                                    ...glassButton,
                                    width: '100%',
                                    display: 'flex',
                                    alignItems: 'flex-start',
                                    gap: 10,
                                    padding: '10px 14px',
                                    fontSize: 12,
                                    color: currentFile === m ? g.accent : g.textDim,
                                    background: currentFile === m ? g.accentGlow : g.glassBg,
                                    cursor: 'pointer',
                                    textAlign: 'left'
                                }}>
                                    <FileText size={14} />
                                    <span
                                        title={m}
                                        style={{
                                            flex: 1,
                                            whiteSpace: 'normal',
                                            overflowWrap: 'anywhere',
                                            wordBreak: 'break-word',
                                            lineHeight: 1.25
                                        }}
                                    >
                                        {m}
                                    </span>
                                </button>
                            ))}
                        </div>

                        <div style={{ marginTop: 12, display: 'flex', flexDirection: 'column', gap: 8 }}>
                            <h3 style={{ fontSize: 10, fontWeight: 700, color: g.textMuted, letterSpacing: '0.15em', textTransform: 'uppercase', margin: 0 }}>Sample File</h3>
                            <input
                                ref={sampleFileInputRef}
                                type="file"
                                accept=".fcs,.FCS"
                                style={{ display: 'none' }}
                                onChange={e => void handleSampleFilePicked(e)}
                            />
                            <select
                                value={currentSample}
                                onChange={e => void handleSampleChange(e.target.value)}
                                style={{
                                    width: '100%',
                                    background: g.inputBg,
                                    color: g.text,
                                    border: `1px solid ${g.glassBorder}`,
                                    borderRadius: 8,
                                    padding: '8px 10px',
                                    fontSize: 12
                                }}
                            >
                                {sampleFiles.length === 0 && <option value="">(no samples found)</option>}
                                {sampleFiles.map(s => (
                                    <option key={s} value={s}>{s}</option>
                                ))}
                            </select>
                            <button
                                onClick={triggerSampleFilePicker}
                                disabled={sampleImportStatus === 'importing'}
                                style={{
                                    ...glassButton,
                                    padding: '8px 10px',
                                    fontSize: 12,
                                    color: g.text,
                                    cursor: sampleImportStatus === 'importing' ? 'default' : 'pointer',
                                    opacity: sampleImportStatus === 'importing' ? 0.7 : 1
                                }}
                            >
                                {sampleImportStatus === 'importing' ? 'Importing...' : 'Open File...'}
                            </button>
                            {sampleImportStatus !== 'idle' && (
                                <span
                                    style={{
                                        fontSize: 11,
                                        color: sampleImportStatus === 'error' ? '#ef4444' : g.textMuted
                                    }}
                                >
                                    {sampleImportMessage}
                                </span>
                            )}
                        </div>

                        <div style={{ marginTop: 12, display: 'flex', flexDirection: 'column', gap: 8 }}>
                            <h3 style={{ fontSize: 10, fontWeight: 700, color: g.textMuted, letterSpacing: '0.15em', textTransform: 'uppercase', margin: 0 }}>View Config</h3>
                            <input
                                value={currentConfig}
                                onChange={e => setCurrentConfig(e.target.value)}
                                placeholder="gui_config.json"
                                style={{
                                    width: '100%',
                                    background: g.inputBg,
                                    color: g.text,
                                    border: `1px solid ${g.glassBorder}`,
                                    borderRadius: 8,
                                    padding: '8px 10px',
                                    fontSize: 12
                                }}
                            />
                            <select
                                value={configFiles.includes(currentConfig) ? currentConfig : ''}
                                onChange={e => setCurrentConfig(e.target.value)}
                                style={{
                                    width: '100%',
                                    background: g.inputBg,
                                    color: g.text,
                                    border: `1px solid ${g.glassBorder}`,
                                    borderRadius: 8,
                                    padding: '8px 10px',
                                    fontSize: 12
                                }}
                            >
                                <option value="">Select saved config...</option>
                                {configFiles.map(c => (
                                    <option key={c} value={c}>{c}</option>
                                ))}
                            </select>
                            <div style={{ display: 'flex', gap: 8 }}>
                                <button
                                    onClick={() => void loadViewConfig(currentConfig)}
                                    style={{
                                        ...glassButton,
                                        flex: 1,
                                        padding: '8px 10px',
                                        fontSize: 12,
                                        color: g.text,
                                        cursor: 'pointer'
                                    }}
                                >
                                    Load
                                </button>
                                <button
                                    onClick={() => void saveViewConfig()}
                                    style={{
                                        ...glassButton,
                                        flex: 1,
                                        padding: '8px 10px',
                                        fontSize: 12,
                                        color: g.text,
                                        cursor: 'pointer'
                                    }}
                                >
                                    Save
                                </button>
                            </div>
                            <span style={{ fontSize: 11, color: g.textMuted }}>
                                {configStatus === 'saving' ? 'Saving config...' :
                                    configStatus === 'loading' ? 'Loading config...' :
                                        configStatus === 'saved' ? 'Config saved.' : ''}
                            </span>
                        </div>
                    </div>
                    <div style={{ flex: 1, overflowY: 'auto', padding: '0 16px 16px' }}>
                        <h3 style={{ fontSize: 10, fontWeight: 700, color: g.textMuted, letterSpacing: '0.15em', textTransform: 'uppercase', marginBottom: 10 }}>Signatures</h3>
                        <div style={{ display: 'flex', flexDirection: 'column', gap: 2 }}>
                            {matrix.map((m, idx) => (
                                <label key={m.Marker} style={{
                                    display: 'flex',
                                    alignItems: 'center',
                                    justifyContent: 'space-between',
                                    padding: '8px 14px',
                                    borderRadius: 10,
                                    cursor: 'pointer',
                                    fontSize: 12,
                                    background: selectedMarkers.includes(m.Marker) ? g.glassBg : 'transparent',
                                    color: selectedMarkers.includes(m.Marker) ? g.text : g.textMuted,
                                    transition: 'all 0.15s ease'
                                }}>
                                    <div style={{ display: 'flex', alignItems: 'center', gap: 10 }}>
                                        <input type="checkbox" style={{ display: 'none' }} checked={selectedMarkers.includes(m.Marker)} onChange={() => {
                                            setSelectedMarkers(prev => prev.includes(m.Marker) ? prev.filter(x => x !== m.Marker) : [...prev, m.Marker]);
                                        }} />
                                        <div style={{
                                            width: 10,
                                            height: 10,
                                            borderRadius: '50%',
                                            backgroundColor: colors[idx % colors.length],
                                            opacity: selectedMarkers.includes(m.Marker) ? 1 : 0.3,
                                            boxShadow: selectedMarkers.includes(m.Marker) ? `0 0 8px ${colors[idx % colors.length]}` : 'none'
                                        }} />
                                        <span>{m.Marker}</span>
                                    </div>
                                    {selectedMarkers.includes(m.Marker) && <Check size={12} color={g.accent} />}
                                </label>
                            ))}
                        </div>
                    </div>
                </aside>

                {/* Main Content */}
                <div style={{ flex: 1, display: 'flex', flexDirection: 'column', padding: 16, gap: 16, overflow: pageScroll ? 'visible' : 'auto' }}>
                    {/* Controls */}
                    {showControls && (
                        <div style={{ ...glassCard, padding: 16, flexShrink: 0 }}>
                            <div style={{ display: 'flex', flexWrap: 'wrap', gap: '16px 32px', fontSize: 12 }}>
                                <div style={{ display: 'flex', alignItems: 'center', gap: 12 }}>
                                    <span style={{ fontWeight: 700, color: g.textMuted, textTransform: 'uppercase', letterSpacing: '0.1em' }}>Sig:</span>
                                    <span style={{ color: g.textMuted }}>W</span>
                                    <input type="range" min="0.2" max="2" step="0.1" value={lineWidth} onChange={e => setLineWidth(Number(e.target.value))} style={{ width: 50 }} />
                                    <span style={{ color: g.textMuted, width: 20 }}>{lineWidth}</span>
                                    <span style={{ color: g.textMuted }}>Op</span>
                                    <input type="range" min="0.1" max="1" step="0.05" value={lineOpacity} onChange={e => setLineOpacity(Number(e.target.value))} style={{ width: 50 }} />
                                    <span style={{ color: g.textMuted, width: 30 }}>{Math.round(lineOpacity * 100)}%</span>
                                    <Palette size={14} color={g.textMuted} />
                                    <select value={colorPalette} onChange={e => setColorPalette(e.target.value as keyof typeof COLOR_PALETTES)} style={{
                                        background: g.inputBg,
                                        color: g.text,
                                        padding: '6px 12px',
                                        borderRadius: 8,
                                        border: `1px solid ${g.glassBorder}`,
                                        fontSize: 12
                                    }}>
                                        <option value="default">Default</option>
                                        <option value="neon">Neon</option>
                                        <option value="pastel">Pastel</option>
                                        <option value="warm">Warm</option>
                                        <option value="cool">Cool</option>
                                    </select>
                                    <span style={{ color: g.textMuted }}>H</span>
                                    <input type="range" min="200" max="500" step="25" value={signatureHeight} onChange={e => setSignatureHeight(Number(e.target.value))} style={{ width: 50 }} />
                                    <span style={{ color: g.textMuted, width: 30 }}>{signatureHeight}</span>
                                    <span style={{ color: g.textMuted }}>DetW</span>
                                    <input type="range" min="10" max="30" step="2" value={signatureDetWidth} onChange={e => setSignatureDetWidth(Number(e.target.value))} style={{ width: 40 }} />
                                    <span style={{ color: g.textMuted, width: 20 }}>{signatureDetWidth}</span>
                                </div>
                                <div style={{ display: 'flex', alignItems: 'center', gap: 12 }}>
                                    <span style={{ fontWeight: 700, color: g.textMuted, textTransform: 'uppercase', letterSpacing: '0.1em' }}>Res:</span>
                                    <span style={{ color: g.textMuted }}>Pt</span>
                                    <input type="range" min="0.5" max="4" step="0.25" value={pointSize} onChange={e => setPointSize(Number(e.target.value))} style={{ width: 40 }} />
                                    <span style={{ color: g.textMuted, width: 20 }}>{pointSize}</span>
                                    <span style={{ color: g.textMuted }}>Op</span>
                                    <input type="range" min="0.1" max="1" step="0.05" value={pointOpacity} onChange={e => setPointOpacity(Number(e.target.value))} style={{ width: 40 }} />
                                    <span style={{ color: g.textMuted, width: 30 }}>{Math.round(pointOpacity * 100)}%</span>
                                    <Circle size={12} color={g.textMuted} />
                                    <div style={{ display: 'flex', gap: 4 }}>
                                        {SCATTER_COLORS.map(c => (
                                            <button key={c} onClick={() => setPointColor(c)} style={{
                                                width: 18,
                                                height: 18,
                                                borderRadius: '50%',
                                                backgroundColor: c,
                                                border: pointColor === c ? '2px solid white' : '2px solid transparent',
                                                cursor: 'pointer',
                                                boxShadow: pointColor === c ? `0 0 8px ${c}` : 'none'
                                            }} />
                                        ))}
                                    </div>
                                    <span style={{ color: g.textMuted }}>Cell</span>
                                    <input type="range" min="80" max="200" step="10" value={residualCellSize} onChange={e => setResidualCellSize(Number(e.target.value))} style={{ width: 40 }} />
                                    <span style={{ color: g.textMuted, width: 30 }}>{residualCellSize}</span>
                                </div>
                            </div>
                        </div>
                    )}

                    {/* Signature Editor */}
                    <div style={{ flexShrink: 0 }}>
                        <div style={{ display: 'flex', alignItems: 'center', gap: 10, marginBottom: 8 }}>
                            <div style={{
                                padding: 8,
                                background: 'rgba(99, 102, 241, 0.15)',
                                borderRadius: 10
                            }}>
                                <Settings2 size={16} color="#818cf8" />
                            </div>
                            <h2 style={{ fontWeight: 700, fontSize: 15, margin: 0 }}>Signature Editor</h2>
                            <span style={{ fontSize: 12, color: g.textMuted }}>({selectedMarkers.length} / {detectors.length})</span>
                        </div>
                        <div style={{ ...glassCard, overflow: 'auto', maxHeight: signatureHeight + 80 }}>
                            <div style={{ width: chartWidth }}>
                                <svg
                                    ref={svgRef}
                                    style={{ cursor: 'crosshair', display: 'block' }}
                                    width={chartWidth}
                                    height={signatureHeight}
                                    viewBox={`0 0 ${chartWidth} 100`}
                                    preserveAspectRatio="none"
                                    onMouseDown={() => setIsDragging(true)}
                                    onMouseUp={() => setIsDragging(false)}
                                    onMouseLeave={() => setIsDragging(false)}
                                    onMouseMove={handleMouseMove}
                                >
                                    {[0, 0.25, 0.5, 0.75, 1].map(tick => (
                                        <line key={tick} x1="0" y1={5 + (1 - tick) * 90} x2={chartWidth} y2={5 + (1 - tick) * 90} stroke={g.gridLine} strokeWidth="0.3" strokeDasharray="3 3" />
                                    ))}
                                    {selectedMarkers.map((m, mIdx) => {
                                        const row = matrix.find(r => r.Marker === m);
                                        if (!row) return null;
                                        const points = detectors.map((det, dIdx) => {
                                            const px = dIdx * signatureDetWidth + signatureDetWidth / 2;
                                            const py = 5 + (1 - Number(row[det])) * 90;
                                            return `${px},${py}`;
                                        }).join(' ');
                                        return (
                                            <polyline key={m} points={points} fill="none" stroke={colors[mIdx % colors.length]} strokeWidth={lineWidth} strokeLinejoin="round" strokeLinecap="round" style={{ opacity: lineOpacity }} />
                                        );
                                    })}
                                </svg>
                                <div style={{
                                    display: 'flex',
                                    borderTop: `1px solid ${g.glassBorder}`,
                                    background: g.inputBg,
                                    width: chartWidth
                                }}>
                                    {detectorLabels.map((label, idx) => (
                                        <div key={idx} style={{ position: 'relative', width: signatureDetWidth, height: 60 }}>
                                            <span style={{
                                                position: 'absolute',
                                                fontSize: 8,
                                                color: g.textMuted,
                                                whiteSpace: 'nowrap',
                                                left: signatureDetWidth / 2,
                                                top: 4,
                                                transform: 'rotate(90deg) translateX(0) translateY(-50%)',
                                                transformOrigin: 'top left'
                                            }}>
                                                {label}
                                            </span>
                                        </div>
                                    ))}
                                </div>
                            </div>
                        </div>
                    </div>

                    {/* Residual Monitor */}
                    <div style={{ flex: pageScroll ? 'none' : 1, minHeight: 0 }}>
                        <div style={{ display: 'flex', alignItems: 'center', gap: 10, marginBottom: 8, flexWrap: 'wrap' }}>
                            <div style={{
                                padding: 8,
                                background: 'rgba(16, 185, 129, 0.15)',
                                borderRadius: 10
                            }}>
                                <Activity size={16} color="#34d399" />
                            </div>
                            <h2 style={{ fontWeight: 700, fontSize: 15, margin: 0 }}>Residual Monitor</h2>
                            <span style={{ fontSize: 12, color: g.textMuted }}>({lowerTriangleCells} cells, {unmixedData.length} events)</span>
                            <span style={{ fontSize: 11, color: '#fb7185', fontWeight: 600, marginLeft: 16, display: 'flex', alignItems: 'center', gap: 6 }}>
                                <Sliders size={12} /> Drag on plots to adjust crosstalk
                            </span>
                            <span style={{ fontSize: 12, color: g.textMuted, marginLeft: 8 }}>Drag Sensitivity:</span>
                            <input
                                type="range"
                                min="0.01"
                                max="0.1"
                                step="0.01"
                                value={dragSensitivity}
                                onChange={e => setDragSensitivity(Number(e.target.value))}
                                style={{ width: 110 }}
                            />
                            <span style={{ width: 42, fontSize: 12, color: g.text, fontVariantNumeric: 'tabular-nums' }}>
                                {dragSensitivity.toFixed(2)}
                            </span>
                        </div>
                        <div style={{ ...glassCard, padding: 12, overflow: 'auto', maxHeight: pageScroll ? 'none' : 'calc(100% - 40px)' }}>
                            <div style={{ display: 'inline-block' }}>
                                {/* Header row for first x-axis label */}
                                <div style={{ display: selectedMarkers.includes(markerNames[0]) ? 'flex' : 'none', alignItems: 'center' }}>
                                    <div style={{ width: 60 }}></div>
                                    <div style={{
                                        width: residualCellSize,
                                        height: 24,
                                        display: 'flex',
                                        alignItems: 'center',
                                        justifyContent: 'center',
                                        fontSize: 10,
                                        color: g.textMuted,
                                        fontWeight: 500
                                    }}>
                                        {markerNames[0]}
                                    </div>
                                </div>
                                {/* Lower triangle rows */}
                                {markerNames.map((rowName, rowIdx) => {
                                    if (rowIdx === 0) return null;
                                    const rowSelected = selectedMarkers.includes(rowName);
                                    const isLastRow = rowIdx === markerNames.length - 1;
                                    return (
                                        <div key={`row-${rowIdx}`} style={{ display: rowSelected ? 'flex' : 'none', alignItems: 'center' }}>
                                            <div style={{
                                                width: 60,
                                                fontSize: 10,
                                                color: g.textMuted,
                                                textAlign: 'right',
                                                fontWeight: 500,
                                                paddingRight: 8,
                                                overflow: 'hidden',
                                                textOverflow: 'ellipsis',
                                                whiteSpace: 'nowrap'
                                            }}>{rowName}</div>
                                            {markerNames.slice(0, rowIdx).map((colName, colIdx) => {
                                                const colSelected = selectedMarkers.includes(colName);
                                                return (
                                                    <div key={`cell-${rowIdx}-${colIdx}`} style={{
                                                        width: residualCellSize,
                                                        height: residualCellSize,
                                                        display: colSelected ? 'block' : 'none',
                                                        border: `1px solid ${g.glassBorder}`,
                                                        background: g.glassBg,
                                                        borderRadius: 8,
                                                        padding: 3,
                                                        margin: 1
                                                    }}>
                                                        <ResidualPlot
                                                            xKey={colName}
                                                            yKey={rowName}
                                                            data={unmixedData}
                                                            size={residualCellSize - 10}
                                                            pointColor={pointColor}
                                                            pointOpacity={pointOpacity}
                                                            pointSize={pointSize}
                                                            sensitivity={dragSensitivity}
                                                            onAdjust={handleResidualAdjust}
                                                        />
                                                    </div>
                                                );
                                            })}
                                            {!isLastRow && (
                                                <div style={{
                                                    width: residualCellSize,
                                                    height: residualCellSize,
                                                    display: 'flex',
                                                    alignItems: 'center',
                                                    justifyContent: 'center',
                                                    fontSize: 10,
                                                    color: g.textMuted,
                                                    fontWeight: 500
                                                }}>
                                                    {markerNames[rowIdx]}
                                                </div>
                                            )}
                                        </div>
                                    );
                                })}
                            </div>
                        </div>
                    </div>
                </div>
            </div>
        </div>
    );
};

export default App;
