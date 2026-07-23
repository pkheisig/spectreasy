import { useEffect, useRef, useState } from 'react';
import { Info, LoaderCircle, Moon, Save, Settings, Sun, X } from 'lucide-react';
import axios from 'axios';
import { ModuleLoadingState } from './ModuleLoadingState';
import ResidualPlot from './ResidualPlot';
import { DetectorSpectrumAxis } from './DetectorSpectrumAxis';
import { DETECTOR_AXIS_FOOTER_HEIGHT, detectorAxisChartWidth } from './detectorAxis';
import { appletCacheKey, loadCachedAppletData, setCachedAppletData } from './appletDataCache';
import {
    API_BASE,
    SIGNATURE_COLORS,
    StyledDropdown,
    UNMIXING_METHODS,
    alignDetectorLabels,
    asScalarString,
    asUnmixingMethod,
    isUnmixingFilename,
    pickPreferredMatrix,
    projectBody,
    projectUrl,
} from './matrixAdjustmentShared';
import type {
    DataRow,
    MatrixPayload,
    MatrixRow,
    UnmixingMethod,
    ViewConfig,
} from './matrixAdjustmentShared';

type MatrixAdjustmentProps = {
    embedded?: boolean;
    cockpitTheme?: 'light' | 'dark' | null;
    projectPath?: string;
    projectRevision?: string;
    initialMatrixFiles?: string[];
    initialSampleFiles?: string[];
    initialUnmixingMethod?: string;
    onRequestExit?: () => void;
};

const App = ({ embedded = false, cockpitTheme = null, projectPath = '', projectRevision = 'empty', initialMatrixFiles, initialSampleFiles, initialUnmixingMethod, onRequestExit }: MatrixAdjustmentProps) => {
    const initialMethod = asUnmixingMethod(initialUnmixingMethod);
    const [matrices, setMatrices] = useState<string[]>(() => initialMatrixFiles ?? []);
    const [currentFile, setCurrentFile] = useState('');
    const [sampleFiles, setSampleFiles] = useState<string[]>(() => initialSampleFiles ?? []);
    const [currentSample, setCurrentSample] = useState('');
    const [matrix, setMatrix] = useState<MatrixRow[]>([]);
    const [loading, setLoading] = useState(true);
    const [previewLoading, setPreviewLoading] = useState(false);
    const [detectors, setDetectors] = useState<string[]>([]);
    const [detectorLabels, setDetectorLabels] = useState<string[]>([]);
    const [unmixedData, setUnmixedData] = useState<DataRow[]>([]);
    const [rawData, setRawData] = useState<DataRow[]>([]);
    const [isUnmixingMatrix, setIsUnmixingMatrix] = useState(false);
    const [saveStatus, setSaveStatus] = useState<'idle' | 'saving' | 'saved'>('idle');
    const [errorMessage, setErrorMessage] = useState('');
    const [guiStateLoaded, setGuiStateLoaded] = useState(false);
    const [unmixingMethod, setUnmixingMethod] = useState<UnmixingMethod>(initialMethod);

    const [residualCellSize, setResidualCellSize] = useState(130);
    const [pointSize, setPointSize] = useState(1.5);
    const [pointOpacity, setPointOpacity] = useState(0.5);
    const [dragSensitivity, setDragSensitivity] = useState(0.1);

    const [theme, setTheme] = useState<'dark' | 'light'>(() => {
        if (embedded && cockpitTheme) return cockpitTheme;
        const stored = localStorage.getItem('spectreasy-theme');
        if (stored === 'dark' || stored === 'light') return stored;
        return window.matchMedia('(prefers-color-scheme: dark)').matches ? 'dark' : 'light';
    });
    const [settingsOpen, setSettingsOpen] = useState(false);
    const [hoveredSignature, setHoveredSignature] = useState<string | null>(null);
    const unmixingMethodRef = useRef<UnmixingMethod>(initialMethod);
    const bootPromiseRef = useRef<Promise<void> | null>(null);
    const unmixRequestIdRef = useRef(0);

    useEffect(() => {
        const previousTitle = document.title;
        document.title = 'Spectreasy · Matrix Adjustment';
        return () => { document.title = previousTitle; };
    }, []);

    const fetchMatrices = async (force = false) => {
        if (!force && initialMatrixFiles !== undefined) {
            setMatrices(initialMatrixFiles);
            return initialMatrixFiles;
        }
        const res = await axios.get(projectUrl('/matrices'));
        const list = Array.isArray(res.data) ? res.data : [];
        setMatrices(list);
        return list;
    };

    const fetchSamples = async () => {
        if (initialSampleFiles !== undefined) {
            setSampleFiles(initialSampleFiles);
            return initialSampleFiles;
        }
        const res = await axios.get(projectUrl('/samples'));
        const list = Array.isArray(res.data) ? res.data : [];
        setSampleFiles(list);
        return list;
    };

    const fetchStatus = async () => {
        if (initialUnmixingMethod) {
            const method = asUnmixingMethod(initialUnmixingMethod);
            unmixingMethodRef.current = method;
            setUnmixingMethod(method);
            return method;
        }
        const result = await axios.get(`${API_BASE}/status`).catch(() => null);
        const method = asUnmixingMethod(result?.data?.unmixing_method);
        unmixingMethodRef.current = method;
        setUnmixingMethod(method);
        return method;
    };

    const fetchUserGuiState = async () => {
        const result = await axios.get(projectUrl('/gui_state?module=matrix_tuner')).catch(() => null);
        if (result?.data?.config) {
            applyConfig(result.data.config);
        }
        setGuiStateLoaded(true);
    };

    const fetchSampleData = async (sampleName: string, matrixData: MatrixRow[], filenameForType: string, detNames: string[]) => {
        const q = sampleName && sampleName.length > 0
            ? `?sample_name=${encodeURIComponent(sampleName)}`
            : '';
        const sampleKey = appletCacheKey('matrix-sample', projectPath, projectRevision, sampleName || '(default)');
        const sampleData = await loadCachedAppletData(sampleKey, async () => {
            const response = await axios.get(projectUrl(`/data${q}`));
            return response.data as Record<string, unknown>;
        });
        if (sampleData.error) {
            setRawData([]);
            setUnmixedData([]);
            return;
        }
        const raw = Array.isArray(sampleData.raw_data) ? sampleData.raw_data as DataRow[] : [];
        setRawData(raw);
        setDetectorLabels(alignDetectorLabels(detNames, sampleData));
        if (sampleData.sample_name) setCurrentSample(asScalarString(sampleData.sample_name));
        void runUnmix(matrixData, raw, filenameForType);
    };

    const fetchData = async (filename = currentFile, sampleName = currentSample) => {
        if (!filename) {
            setLoading(false);
            return;
        }
        setLoading(true);
        setErrorMessage('');
        try {
            const matrixKey = appletCacheKey('matrix-file', projectPath, projectRevision, filename);
            const matrixData = await loadCachedAppletData(matrixKey, async () => {
                const response = await axios.get(projectUrl(`/load_matrix?filename=${encodeURIComponent(filename)}`));
                if (response.data?.error) throw new Error(asScalarString(response.data.error));
                return Array.isArray(response.data) ? response.data as MatrixRow[] : [];
            });
            if (matrixData.length === 0) {
                throw new Error(`Matrix ${filename} contains no rows.`);
            }
            setMatrix(matrixData);
            setIsUnmixingMatrix(isUnmixingFilename(filename));
            const detNames = Object.keys(matrixData[0]).filter(k => k !== 'Marker');
            setDetectors(detNames);
            setDetectorLabels(detNames);
            const activeSample = sampleName && sampleName.length > 0 ? sampleName : (sampleFiles[0] || '');
            await fetchSampleData(activeSample, matrixData, filename, detNames);
            setCurrentFile(filename);
        } catch (error) {
            setMatrix([]);
            setDetectors([]);
            setDetectorLabels([]);
            setUnmixedData([]);
            setErrorMessage(error instanceof Error ? error.message : 'Could not load the selected matrix.');
        } finally {
            setLoading(false);
        }
    };

    useEffect(() => {
        const init = async () => {
            try {
                const [, mats, samples] = await Promise.all([
                    fetchStatus(),
                    fetchMatrices(),
                    fetchSamples(),
                    fetchUserGuiState(),
                ]);
                const firstMatrix = pickPreferredMatrix(mats, currentFile);
                const firstSample = samples[0] || '';
                if (firstSample) setCurrentSample(firstSample);
                if (firstMatrix) await fetchData(firstMatrix, firstSample);
                else setLoading(false);
            } catch (error) {
                setErrorMessage(error instanceof Error ? error.message : 'Matrix Adjustment could not reach the local R backend.');
                setLoading(false);
            }
        };
        if (!bootPromiseRef.current) bootPromiseRef.current = init();
        void bootPromiseRef.current;
        // Initial boot should run once; the called helpers intentionally read their current defaults.
        // eslint-disable-next-line react-hooks/exhaustive-deps
    }, []);

    async function runUnmix(currentM: MatrixRow[], currentRaw: DataRow[], filename: string | boolean = currentFile) {
        const requestId = ++unmixRequestIdRef.current;
        if (!Array.isArray(currentM) || currentM.length === 0 || !Array.isArray(currentRaw) || currentRaw.length === 0) {
            setUnmixedData([]);
            setPreviewLoading(false);
            return false;
        }
        setPreviewLoading(true);
        const M_obj: MatrixPayload = {};
        currentM.forEach(row => {
            M_obj[row.Marker] = { ...row };
            delete M_obj[row.Marker].Marker;
        });
        let useType = 'reference';
        if (typeof filename === 'string') {
            if (isUnmixingFilename(filename)) useType = 'unmixing';
        } else if (typeof filename === 'boolean') {
            useType = filename ? 'unmixing' : 'reference';
        }
        try {
            const res = await axios.post(`${API_BASE}/unmix`, projectBody({
                matrix_json: M_obj,
                raw_data_json: currentRaw,
                type: useType,
                matrix_filename: typeof filename === 'string' ? filename : currentFile,
                method: unmixingMethodRef.current
            }));
            if (requestId !== unmixRequestIdRef.current) return false;
            if (Array.isArray(res.data)) {
                setUnmixedData(res.data as DataRow[]);
                return true;
            }
            throw new Error(asScalarString(res.data?.error, 'The R backend returned no unmixed data.'));
        } catch (error) {
            if (requestId === unmixRequestIdRef.current) {
                setUnmixedData([]);
                setErrorMessage(error instanceof Error ? error.message : 'Sample unmixing failed.');
            }
        } finally {
            if (requestId === unmixRequestIdRef.current) setPreviewLoading(false);
        }
        return false;
    }

    const handleUnmixingMethodChange = (method: UnmixingMethod) => {
        unmixingMethodRef.current = method;
        setUnmixingMethod(method);
        if (matrix.length > 0 && rawData.length > 0) {
            void runUnmix(matrix, rawData, isUnmixingMatrix);
        }
    };

    const sampleHasDetectorChannels = (sampleRows = rawData) => {
        if (!Array.isArray(sampleRows) || sampleRows.length === 0 || detectors.length === 0) return false;
        const sampleColumns = new Set(Object.keys(sampleRows[0]));
        return detectors.some(detector => sampleColumns.has(detector));
    };

    const handleResidualAdjust = (xMarker: string, yMarker: string, alpha: number) => {
        if (!sampleHasDetectorChannels()) return;

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
        runUnmix(newMatrix, rawData, isUnmixingMatrix);
    };

    const saveMatrix = async () => {
        if (!currentFile || matrix.length === 0) return;
        setSaveStatus('saving');
        setErrorMessage('');
        const newName = /\.csv$/i.test(currentFile)
            ? currentFile.replace(/\.csv$/i, '_adjusted.csv')
            : `${currentFile}_adjusted.csv`;
        try {
            const result = await axios.post(`${API_BASE}/save_matrix`, projectBody({
                filename: newName,
                source_filename: currentFile,
                matrix_json: matrix
            }));
            if (result.data?.error) throw new Error(asScalarString(result.data.error));
            setSaveStatus('saved');
            const savedMatrixKey = appletCacheKey('matrix-file', projectPath, projectRevision, newName);
            setCachedAppletData(savedMatrixKey, matrix);
            await fetchMatrices(true);
            window.setTimeout(() => setSaveStatus('idle'), 2000);
        } catch (error) {
            setErrorMessage(error instanceof Error ? error.message : 'The adjusted matrix could not be saved.');
            setSaveStatus('idle');
        }
    };

    function applyConfig(cfg: unknown) {
        if (!cfg || typeof cfg !== 'object') return;
        const viewConfig = cfg as ViewConfig;
        if (typeof viewConfig.residualCellSize === 'number') setResidualCellSize(viewConfig.residualCellSize);
        if (typeof viewConfig.pointSize === 'number') setPointSize(viewConfig.pointSize);
        if (typeof viewConfig.pointOpacity === 'number') setPointOpacity(viewConfig.pointOpacity);
        if (typeof viewConfig.dragSensitivity === 'number') {
            setDragSensitivity(Math.max(0, Math.min(0.3, viewConfig.dragSensitivity)));
        }
        if (!embedded && (viewConfig.theme === 'dark' || viewConfig.theme === 'light')) setTheme(viewConfig.theme);
    }

    const buildConfig = () => ({
        residualCellSize,
        pointSize,
        pointOpacity,
        dragSensitivity,
        ...(!embedded ? { theme } : {})
    });

    useEffect(() => {
        if (!guiStateLoaded) return;
        const timer = window.setTimeout(() => {
            void axios.post(`${API_BASE}/gui_state`, projectBody({
                module: 'matrix_tuner',
                config_json: buildConfig()
            })).catch(() => null);
        }, 500);
        return () => window.clearTimeout(timer);
        // eslint-disable-next-line react-hooks/exhaustive-deps
    }, [residualCellSize, pointSize, pointOpacity, dragSensitivity, theme, guiStateLoaded]);

    const handleSampleChange = async (sampleName: string) => {
        setCurrentSample(sampleName);
        if (!sampleName || matrix.length === 0) return;
        setLoading(true);
        setErrorMessage('');
        try {
            await fetchSampleData(sampleName, matrix, currentFile, detectors);
        } catch (error) {
            setErrorMessage(error instanceof Error ? error.message : 'Could not load the selected sample.');
        } finally {
            setLoading(false);
        }
    };

    const svgRef = useRef<SVGSVGElement>(null);
    const colors = SIGNATURE_COLORS;
    const markerNames = matrix.map(m => m.Marker);
    const selectedMarkers = markerNames;
    const signatureDetectorLabels = detectorLabels.length === detectors.length ? detectorLabels : detectors;
    const chartWidth = detectorAxisChartWidth(detectors.length);
    const chartHeight = 230;
    const spectrumLeft = 42;
    const spectrumRight = chartWidth - 8;
    const spectrumPlotWidth = spectrumRight - spectrumLeft;
    const detectorX = (idx: number) => {
        if (detectors.length <= 1) return spectrumLeft + spectrumPlotWidth / 2;
        return spectrumLeft + (idx / (detectors.length - 1)) * spectrumPlotWidth;
    };
    const signaturePath = (row: MatrixRow) => detectors.map((det, idx) => {
        const x = detectorX(idx);
        const y = chartHeight - Number(row[det]) * (chartHeight - 32) - 24;
        return `${idx === 0 ? 'M' : 'L'}${x.toFixed(1)},${y.toFixed(1)}`;
    }).join(' ');
    const detectorAxisEntries = detectors.map((detector, index) => ({
        detector,
        label: signatureDetectorLabels[index] || detector,
    }));
    const detectorGradientId = 'matrix-detector-spectrum';
    const renderedTheme = embedded && cockpitTheme ? cockpitTheme : theme;
    const residualRenderKey = [
        residualCellSize,
        pointSize,
        pointOpacity,
        renderedTheme === 'dark' ? '#ff7a1a' : '#8052c7',
        dragSensitivity,
        unmixedData.length
    ].join(':');
    const canAdjustResiduals = sampleHasDetectorChannels();
    const residualAdjustmentText = canAdjustResiduals
        ? 'Drag on plots to adjust crosstalk'
        : 'Load raw detector FCS to adjust crosstalk';

    // iOS 26 Glassy Theme
    const glassyTheme = renderedTheme === 'dark' ? {
        bgGradient: '#121816',
        // Glass panel styles
        glassBg: '#1e2522',
        glassBorder: '#3b4541',
        glassHighlight: '#252e2a',
        glassBlur: 'none',
        // Text colors
        text: '#f0f3f2',
        textMuted: '#9ba7a2',
        textDim: '#a9b7b1',
        // Accent colors - softer teal/cyan
        accent: '#55c7b6',
        accentGlow: 'rgba(85, 199, 182, 0.15)',
        // Input styles
        inputBg: '#121816',
        // Grid
        gridLine: 'rgba(169, 183, 177, 0.15)',
    } : {
        bgGradient: '#f5f3ee',
        glassBg: '#fffdfa',
        glassBorder: '#d8d5cc',
        glassHighlight: '#f8f7f3',
        glassBlur: 'none',
        text: '#17201d',
        textMuted: '#68736e',
        textDim: '#6d756f',
        accent: '#1f7a6d',
        accentGlow: 'rgba(31, 122, 109, 0.15)',
        inputBg: '#fffdfa',
        gridLine: 'rgba(109, 117, 111, 0.18)',
    };

    const g = glassyTheme;

    useEffect(() => {
        if (embedded) return;
        document.documentElement.style.backgroundColor = g.bgGradient;
        document.body.style.backgroundColor = g.bgGradient;
        document.documentElement.style.setProperty('--bg-app', g.bgGradient);
        document.documentElement.dataset.theme = theme;
        localStorage.setItem('spectreasy-theme', theme);
    }, [embedded, g.bgGradient, theme]);

    const glassCard = {
        background: g.glassBg,
        border: `1px solid ${g.glassBorder}`,
        borderRadius: '8px',
        boxShadow: 'none',
    };

    const glassButton = {
        background: g.glassBg,
        border: `1px solid ${g.glassBorder}`,
        borderRadius: '7px',
    };

    if (loading) return <ModuleLoadingState label="Loading Matrix Adjustment" theme={renderedTheme} />;

    const lowerTriangleCells = selectedMarkers.length * (selectedMarkers.length - 1) / 2;

    return (
        <div style={{
            minHeight: '100vh',
            height: 'auto',
            width: '100vw',
            maxWidth: '100vw',
            display: 'flex',
            flexDirection: 'column',
            backgroundColor: g.bgGradient,
            backgroundImage: renderedTheme === 'dark'
                ? 'linear-gradient(90deg, rgba(255,255,255,.02) 1px, transparent 1px), linear-gradient(rgba(255,255,255,.02) 1px, transparent 1px)'
                : 'linear-gradient(90deg, rgba(38,63,115,.035) 1px, transparent 1px), linear-gradient(rgba(38,63,115,.035) 1px, transparent 1px)',
            backgroundSize: '22px 22px',
            color: g.text,
            fontFamily: 'Avenir Next, "Segoe UI", sans-serif',
            overflowX: 'hidden',
            overflowY: 'auto'
        }}>
            {/* Navbar */}
            <header style={{
                minHeight: 76,
                flexShrink: 0,
                borderBottom: `1px solid ${g.glassBorder}`,
                display: 'flex',
                alignItems: 'center',
                flexWrap: 'wrap',
                gap: 16,
                padding: '10px 18px',
                justifyContent: 'space-between',
                background: g.glassBg,
                position: 'sticky',
                top: 0,
                zIndex: 50
            }}>
                <div style={{ minWidth: 230 }}>
                    <h1 style={{ fontWeight: 750, fontSize: 21, margin: 0, letterSpacing: 0 }}>Matrix Adjustment</h1>
                    <div style={{ display: 'flex', alignItems: 'center', flexWrap: 'wrap', gap: '4px 10px', marginTop: 3, color: g.textMuted, fontSize: 11 }}>
                        <span>{detectors.length} detectors · {markerNames.length} signatures · {lowerTriangleCells} residual plots · {unmixedData.length.toLocaleString()} events</span>
                        <span style={{ display: 'inline-flex', alignItems: 'center', gap: 4, color: canAdjustResiduals ? g.accent : g.textMuted, fontWeight: 650 }}>
                            <Info size={12} /> {residualAdjustmentText}
                        </span>
                    </div>
                </div>
                <div style={{ display: 'flex', alignItems: 'flex-end', flexWrap: 'wrap', gap: 8, flex: '1 1 560px', justifyContent: 'flex-end' }}>
                    <StyledDropdown
                        label="Matrix"
                        value={currentFile}
                        options={matrices}
                        emptyLabel="(no matrices found)"
                        width={390}
                        showEnd
                        theme={g}
                        onChange={value => void fetchData(value, currentSample)}
                    />
                    <StyledDropdown
                        label="Sample"
                        value={currentSample}
                        options={sampleFiles}
                        emptyLabel="(no samples found)"
                        width={210}
                        theme={g}
                        onChange={value => void handleSampleChange(value)}
                    />
                    <StyledDropdown
                        label="Method"
                        value={unmixingMethod}
                        options={UNMIXING_METHODS}
                        width={126}
                        theme={g}
                        onChange={value => handleUnmixingMethodChange(asUnmixingMethod(value))}
                    />
                    <div style={{ display: 'flex', alignItems: 'center', gap: 8, position: 'relative' }}>
                    <button aria-label="Plot settings" title="Plot settings" aria-expanded={settingsOpen} onClick={() => setSettingsOpen(open => !open)} style={{ ...glassButton, width: 38, height: 38, display: 'grid', placeItems: 'center', color: settingsOpen ? g.accent : g.textDim, cursor: 'pointer', background: settingsOpen ? g.glassHighlight : g.glassBg }}>
                        <Settings size={18} />
                    </button>
                    <button aria-label={saveStatus === 'saving' ? 'Saving adjusted matrix' : saveStatus === 'saved' ? 'Adjusted matrix saved' : 'Save adjusted matrix'} title={saveStatus === 'saving' ? 'Saving…' : saveStatus === 'saved' ? 'Saved' : 'Save adjusted matrix'} disabled={saveStatus === 'saving'} onClick={saveMatrix} style={{
                        background: g.accent,
                        border: `1px solid ${g.accent}`,
                        borderRadius: 7,
                        width: 38,
                        height: 38,
                        color: 'white',
                        cursor: 'pointer',
                        display: 'grid',
                        placeItems: 'center',
                        boxShadow: 'none'
                    }}>
                        <Save size={17} />
                    </button>
                    {embedded && onRequestExit && <button
                        type="button"
                        onClick={onRequestExit}
                        aria-label="Close matrix adjustment and return to cockpit"
                        autoFocus
                        style={{ ...glassButton, minWidth: 78, height: 38, padding: '0 11px', color: g.text, cursor: 'pointer', display: 'inline-flex', alignItems: 'center', justifyContent: 'center', gap: 6, fontSize: 12, fontWeight: 700 }}
                    >
                        <X size={16} /> Close
                    </button>}
                    {settingsOpen && (
                        <div role="dialog" aria-label="Settings" style={{
                            position: 'absolute',
                            top: 46,
                            right: 46,
                            width: 330,
                            padding: 16,
                            border: `1px solid ${g.glassBorder}`,
                            borderRadius: 9,
                            background: g.glassBg,
                            boxShadow: '0 18px 44px rgba(23, 32, 29, 0.22)',
                            color: g.text,
                            zIndex: 80
                        }}>
                            <div style={{ marginBottom: 14 }}>
                                <strong style={{ display: 'block', fontSize: 14 }}>Settings</strong>
                                <span style={{ display: 'block', marginTop: 3, color: g.textMuted, fontSize: 11 }}>Changes are saved automatically</span>
                            </div>
                            {!embedded && <div style={{ display: 'flex', alignItems: 'center', justifyContent: 'space-between', padding: '9px 0 12px', borderTop: `1px solid ${g.glassBorder}` }}>
                                <span style={{ color: g.textDim, fontSize: 12, fontWeight: 700 }}>Appearance</span>
                                <button aria-label={theme === 'dark' ? 'Use light mode' : 'Use dark mode'} onClick={() => setTheme(theme === 'dark' ? 'light' : 'dark')} style={{ ...glassButton, minWidth: 92, height: 32, padding: '0 10px', color: g.text, cursor: 'pointer', display: 'flex', alignItems: 'center', justifyContent: 'center', gap: 7, fontSize: 11, fontWeight: 700 }}>
                                    {theme === 'dark' ? <><Moon size={14} /> Dark</> : <><Sun size={14} /> Light</>}
                                </button>
                            </div>}
                            {[
                                { label: 'Point size', value: pointSize, display: pointSize.toFixed(2), min: 0.5, max: 4, step: 0.25, set: setPointSize },
                                { label: 'Point opacity', value: pointOpacity, display: `${Math.round(pointOpacity * 100)}%`, min: 0.1, max: 1, step: 0.05, set: setPointOpacity },
                                { label: 'Cell size', value: residualCellSize, display: `${residualCellSize}px`, min: 80, max: 200, step: 10, set: setResidualCellSize },
                                { label: 'Drag sensitivity', value: dragSensitivity, display: dragSensitivity.toFixed(2), min: 0, max: 0.3, step: 0.01, set: setDragSensitivity }
                            ].map(setting => (
                                <label key={setting.label} style={{ display: 'grid', gridTemplateColumns: '104px 1fr 48px', alignItems: 'center', gap: 10, minHeight: 40, borderTop: `1px solid ${g.glassBorder}`, fontSize: 12 }}>
                                    <span style={{ color: g.textDim, fontWeight: 700 }}>{setting.label}</span>
                                    <input type="range" min={setting.min} max={setting.max} step={setting.step} value={setting.value} onChange={event => setting.set(Number(event.target.value))} style={{ width: '100%', accentColor: g.accent }} />
                                    <strong style={{ textAlign: 'right', fontVariantNumeric: 'tabular-nums', fontSize: 11 }}>{setting.display}</strong>
                                </label>
                            ))}
                        </div>
                    )}
                    </div>
                </div>
            </header>
            {errorMessage && (
                <div role="alert" style={{ margin: '12px 18px 0', padding: '10px 12px', color: g.accent, border: `1px solid ${g.accent}`, borderRadius: 7, background: g.glassBg }}>
                    {errorMessage}
                </div>
            )}

            <div style={{ flex: 1, display: 'flex', minWidth: 0, overflowX: 'hidden', overflowY: 'visible' }}>
                <div style={{ flex: 1, minWidth: 0, maxWidth: '100%', display: 'flex', flexDirection: 'column', padding: 16, gap: 16, overflowX: 'hidden', overflowY: 'visible' }}>
                    <div style={{ flexShrink: 0 }}>
                        <div style={{ ...glassCard, overflow: 'hidden', maxWidth: '100%' }}>
                            <div style={{ width: '100%', minWidth: 0, overflowX: 'auto' }}>
                                <svg
                                    ref={svgRef}
                                    style={{ display: 'block', width: chartWidth, minWidth: chartWidth, maxWidth: 'none' }}
                                    width={chartWidth}
                                    height={chartHeight + DETECTOR_AXIS_FOOTER_HEIGHT}
                                    viewBox={`0 0 ${chartWidth} ${chartHeight + DETECTOR_AXIS_FOOTER_HEIGHT}`}
                                    role="img"
                                    aria-label="Reference signatures across detector channels grouped by laser"
                                >
                                    {[0, 25, 50, 75, 100].map(tick => {
                                        const y = chartHeight - (tick / 100) * (chartHeight - 32) - 24;
                                        return (
                                            <g key={tick}>
                                                <line x1={spectrumLeft} y1={y} x2={spectrumRight} y2={y} stroke={g.gridLine} strokeWidth={1} />
                                            </g>
                                        );
                                    })}
                                    <DetectorSpectrumAxis
                                        entries={detectorAxisEntries}
                                        xForIndex={detectorX}
                                        left={spectrumLeft}
                                        right={spectrumRight}
                                        plotTop={6}
                                        baselineY={chartHeight - 24}
                                        gridColor={g.gridLine}
                                        axisColor={g.glassBorder}
                                        textColor={g.textMuted}
                                        gradientId={detectorGradientId}
                                    />
                                    {selectedMarkers.map((m, mIdx) => {
                                        const row = matrix.find(r => r.Marker === m);
                                        if (!row) return null;
                                        const pathData = signaturePath(row);
                                        const color = colors[mIdx % colors.length];
                                        const isHovered = hoveredSignature === m;
                                        const hasHoverActive = hoveredSignature !== null;
                                        return (
                                            <g key={m}>
                                                <path
                                                    d={pathData}
                                                    fill="none"
                                                    stroke={color}
                                                    strokeWidth={isHovered ? 4.2 : 2.4}
                                                    strokeLinejoin="round"
                                                    strokeLinecap="round"
                                                    opacity={hasHoverActive ? (isHovered ? 1.0 : 0.15) : 0.92}
                                                    style={{
                                                        transition: 'all 0.15s ease-in-out',
                                                        pointerEvents: 'none',
                                                        filter: isHovered ? `drop-shadow(0 0 5px ${color})` : 'none'
                                                    }}
                                                />
                                                <path
                                                    d={pathData}
                                                    fill="none"
                                                    stroke="transparent"
                                                    strokeWidth={22}
                                                    strokeLinejoin="round"
                                                    strokeLinecap="round"
                                                    onMouseEnter={() => setHoveredSignature(m)}
                                                    onMouseLeave={() => setHoveredSignature(null)}
                                                    style={{ cursor: 'pointer', pointerEvents: 'stroke' }}
                                                />
                                            </g>
                                        );
                                    })}
                                </svg>
                            </div>
                        </div>
                    </div>

                    <div style={{ flex: 'none', minHeight: 0 }}>
                        <div style={{ ...glassCard, padding: 12, overflow: 'auto', maxHeight: 'none', position: 'relative' }}>
                            {previewLoading && unmixedData.length === 0 ? (
                                <div role="status" aria-live="polite" style={{ minHeight: 220, display: 'flex', alignItems: 'center', justifyContent: 'center', gap: 11, color: g.accent, fontSize: 13, fontWeight: 700 }}>
                                    <LoaderCircle className="module-loading-spinner" size={25} aria-hidden="true" />
                                    Calculating residual preview
                                </div>
                            ) : <>
                            {previewLoading && (
                                <div role="status" aria-live="polite" style={{ position: 'sticky', left: 12, top: 0, zIndex: 5, width: 'fit-content', display: 'flex', alignItems: 'center', gap: 7, marginBottom: 8, padding: '7px 10px', border: `1px solid ${g.glassBorder}`, borderRadius: 7, color: g.accent, background: g.glassBg, fontSize: 11, fontWeight: 700 }}>
                                    <LoaderCircle className="module-loading-spinner" size={16} aria-hidden="true" />
                                    Updating residual preview
                                </div>
                            )}
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
                                                            key={`${rowName}-${colName}-${residualRenderKey}`}
                                                            xKey={colName}
                                                            yKey={rowName}
                                                            data={unmixedData}
                                                            size={residualCellSize - 10}
                                                            pointColor={renderedTheme === 'dark' ? '#ff7a1a' : '#8052c7'}
                                                            pointOpacity={pointOpacity}
                                                            pointSize={pointSize}
                                                            sensitivity={dragSensitivity}
                                                            canAdjust={canAdjustResiduals}
                                                            disabledReason={residualAdjustmentText}
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
                            </>}
                        </div>
                    </div>
                </div>
            </div>
        </div>
    );
};

export default App;
