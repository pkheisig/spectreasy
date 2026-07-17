import { useState, useEffect, useId, useRef } from 'react';
import { Check, ChevronDown, Settings, Save, RefreshCw, Sun, Moon, Info } from 'lucide-react';
import axios from 'axios';
import ResidualPlot from './ResidualPlot';
import { resolveApiBase } from './apiBase';

const API_BASE = resolveApiBase();
const currentProjectPath = () => window.sessionStorage.getItem('spectreasy-project-path') || '';
const projectUrl = (path: string) => {
    const separator = path.includes('?') ? '&' : '?';
    const projectPath = currentProjectPath();
    return `${API_BASE}${path}${projectPath ? `${separator}project_path=${encodeURIComponent(projectPath)}` : ''}`;
};
const projectBody = <T extends Record<string, unknown>>(body: T) => ({ ...body, projectPath: currentProjectPath(), project_path: currentProjectPath() });

interface MatrixRow {
    Marker: string;
    [key: string]: string | number | null | undefined;
}

type DataRow = Record<string, string | number | null | undefined>;
type MatrixPayload = Record<string, Record<string, string | number | null | undefined>>;
const UNMIXING_METHODS = ['Spectreasy', 'AutoSpectral', 'OLS', 'WLS', 'RWLS', 'NNLS'] as const;
type UnmixingMethod = typeof UNMIXING_METHODS[number];
type ViewConfig = {
    residualCellSize?: number;
    pointSize?: number;
    pointOpacity?: number;
    dragSensitivity?: number;
    theme?: 'dark' | 'light';
};

const SIGNATURE_COLORS = ['#60a5fa', '#f472b6', '#34d399', '#fbbf24', '#a78bfa', '#fb7185', '#38bdf8', '#4ade80', '#facc15', '#c084fc'];

const isUnmixingFilename = (filename: string) => {
    const lower = filename.toLowerCase();
    return lower.includes('unmixing') || lower.includes('w_');
};

const pickPreferredMatrix = (files: string[], current: string) => {
    if (files.length === 0) return '';
    if (current && files.includes(current)) return current;

    const preferredNames = [
        'scc_reference_matrix.csv',
        'reference_matrix.csv',
        'refined_reference_matrix.csv',
        'scc_unmixing_matrix.csv',
        'refined_unmixing_matrix.csv'
    ];
    for (const name of preferredNames) {
        const hit = files.find(f => f.toLowerCase() === name);
        if (hit) return hit;
    }

    const referenceHit = files.find(f => f.toLowerCase().includes('reference'));
    if (referenceHit) return referenceHit;

    const unmixingHit = files.find(isUnmixingFilename);
    if (unmixingHit) return unmixingHit;

    return files[0];
};

const asScalarString = (value: unknown, fallback = '') => {
    if (Array.isArray(value)) return asScalarString(value[0], fallback);
    if (typeof value === 'string') return value;
    if (value == null) return fallback;
    return String(value);
};

const asUnmixingMethod = (value: unknown): UnmixingMethod => {
    const text = asScalarString(value, 'Spectreasy');
    return UNMIXING_METHODS.includes(text as UnmixingMethod) ? text as UnmixingMethod : 'Spectreasy';
};

type DropdownTheme = {
    inputBg: string;
    glassBg: string;
    glassBorder: string;
    glassHighlight: string;
    text: string;
    textMuted: string;
    accent: string;
    accentGlow: string;
};

type StyledDropdownProps = {
    label: string;
    value: string;
    options: readonly string[];
    emptyLabel?: string;
    width: number;
    showEnd?: boolean;
    theme: DropdownTheme;
    onChange: (value: string) => void;
};

const StyledDropdown = ({ label, value, options, emptyLabel = '(none available)', width, showEnd = false, theme, onChange }: StyledDropdownProps) => {
    const [open, setOpen] = useState(false);
    const [activeIndex, setActiveIndex] = useState(0);
    const rootRef = useRef<HTMLDivElement>(null);
    const listboxRef = useRef<HTMLDivElement>(null);
    const listboxId = useId();

    useEffect(() => {
        if (!open) return;
        const closeOnOutside = (event: PointerEvent) => {
            if (!rootRef.current?.contains(event.target as Node)) setOpen(false);
        };
        const closeOnEscape = (event: KeyboardEvent) => {
            if (event.key === 'Escape') setOpen(false);
        };
        window.addEventListener('pointerdown', closeOnOutside);
        window.addEventListener('keydown', closeOnEscape);
        return () => {
            window.removeEventListener('pointerdown', closeOnOutside);
            window.removeEventListener('keydown', closeOnEscape);
        };
    }, [open]);

    useEffect(() => {
        if (!open || !showEnd) return;
        const frame = window.requestAnimationFrame(() => {
            if (listboxRef.current) listboxRef.current.scrollLeft = listboxRef.current.scrollWidth;
        });
        return () => window.cancelAnimationFrame(frame);
    }, [open, showEnd]);

    const displayValue = value || emptyLabel;
    const handleKeyDown = (event: React.KeyboardEvent) => {
        if (!options.length) return;
        if (event.key === 'Escape') {
            setOpen(false);
            return;
        }
        if (event.key === 'ArrowDown' || event.key === 'ArrowUp' || event.key === 'Home' || event.key === 'End') {
            event.preventDefault();
            setOpen(true);
            setActiveIndex(previous => {
                if (event.key === 'Home') return 0;
                if (event.key === 'End') return options.length - 1;
                const delta = event.key === 'ArrowDown' ? 1 : -1;
                return (previous + delta + options.length) % options.length;
            });
        } else if ((event.key === 'Enter' || event.key === ' ') && open) {
            event.preventDefault();
            onChange(options[activeIndex]);
            setOpen(false);
        }
    };
    return (
        <div ref={rootRef} onKeyDown={handleKeyDown} style={{ position: 'relative', display: 'flex', flexDirection: 'column', gap: 3, minWidth: 0 }}>
            <span style={{ color: theme.textMuted, fontSize: 11, fontWeight: 700, letterSpacing: '0.1em', textTransform: 'uppercase' }}>{label}</span>
            <button
                type="button"
                aria-haspopup="listbox"
                aria-expanded={open}
                aria-controls={listboxId}
                aria-activedescendant={open && options[activeIndex] ? `${listboxId}-${activeIndex}` : undefined}
                aria-label={`${label}: ${displayValue}`}
                title={displayValue}
                onClick={() => {
                    setActiveIndex(Math.max(0, options.indexOf(value)));
                    setOpen(previous => !previous);
                }}
                style={{
                    width,
                    maxWidth: '32vw',
                    height: 34,
                    minHeight: 34,
                    display: 'flex',
                    alignItems: 'center',
                    justifyContent: 'space-between',
                    gap: 8,
                    padding: '0 9px 0 10px',
                    border: `1px solid ${open ? theme.accent : theme.glassBorder}`,
                    borderRadius: 7,
                    background: open ? theme.glassHighlight : theme.inputBg,
                    color: value ? theme.text : theme.textMuted,
                    boxShadow: open ? `0 0 0 2px ${theme.accentGlow}` : 'none',
                    fontFamily: 'inherit',
                    fontSize: 11,
                    fontWeight: 500,
                    textAlign: 'left',
                    cursor: options.length > 0 ? 'pointer' : 'default'
                }}
            >
                {showEnd ? (
                    <span style={{ minWidth: 0, flex: 1, overflow: 'hidden', display: 'flex', justifyContent: 'flex-end' }}>
                        <span style={{ flex: '0 0 auto', whiteSpace: 'nowrap' }}>{displayValue}</span>
                    </span>
                ) : (
                    <span style={{ minWidth: 0, overflow: 'hidden', textOverflow: 'ellipsis', whiteSpace: 'nowrap' }}>{displayValue}</span>
                )}
                <ChevronDown size={14} style={{ flex: '0 0 auto', transform: open ? 'rotate(180deg)' : 'none', transition: 'transform 140ms ease' }} />
            </button>
            {open && (
                <div
                    ref={listboxRef}
                    id={listboxId}
                    role="listbox"
                    aria-label={label}
                    style={{
                        position: 'absolute',
                        top: 'calc(100% + 6px)',
                        left: showEnd ? 'auto' : 0,
                        right: showEnd ? 0 : 'auto',
                        zIndex: 100,
                        width: showEnd ? Math.max(width, 520) : width,
                        maxWidth: showEnd ? 'min(70vw, 720px)' : 'min(32vw, 360px)',
                        maxHeight: 280,
                        overflow: 'auto',
                        padding: 5,
                        border: `1px solid ${theme.glassBorder}`,
                        borderRadius: 9,
                        background: theme.glassBg,
                        boxShadow: '0 18px 44px rgba(23, 32, 29, 0.24)'
                    }}
                >
                    {options.length === 0 ? (
                        <div style={{ padding: '9px 10px', color: theme.textMuted, fontSize: 11 }}>{emptyLabel}</div>
                    ) : options.map((option, index) => {
                        const selected = option === value;
                        return (
                            <button
                                id={`${listboxId}-${index}`}
                                type="button"
                                role="option"
                                aria-selected={selected}
                                key={option}
                                title={option}
                                onClick={() => {
                                    onChange(option);
                                    setOpen(false);
                                }}
                                style={{
                                    width: showEnd ? 'max-content' : '100%',
                                    minWidth: '100%',
                                    minHeight: 31,
                                    display: 'grid',
                                    gridTemplateColumns: showEnd ? '16px max-content' : '16px minmax(0, 1fr)',
                                    alignItems: 'center',
                                    gap: 7,
                                    padding: '0 8px',
                                    border: 0,
                                    borderRadius: 6,
                                    background: selected || activeIndex === index ? theme.glassHighlight : 'transparent',
                                    color: selected ? theme.accent : theme.text,
                                    fontFamily: 'inherit',
                                    fontSize: 11,
                                    fontWeight: selected ? 600 : 400,
                                    textAlign: 'left',
                                    cursor: 'pointer'
                                }}
                                onMouseEnter={event => { if (!selected) event.currentTarget.style.background = theme.glassHighlight; }}
                                onFocus={() => setActiveIndex(index)}
                                onMouseLeave={event => { if (!selected) event.currentTarget.style.background = 'transparent'; }}
                            >
                                <span style={{ display: 'grid', placeItems: 'center' }}>{selected && <Check size={13} />}</span>
                                <span style={{ minWidth: 0, overflow: showEnd ? 'visible' : 'hidden', textOverflow: showEnd ? 'clip' : 'ellipsis', whiteSpace: 'nowrap' }}>{option}</span>
                            </button>
                        );
                    })}
                </div>
            )}
        </div>
    );
};

const alignDetectorLabels = (detNames: string[], payload: Record<string, unknown>) => {
    const rawNames = Array.isArray(payload.detector_names) ? payload.detector_names.map(v => String(v)) : [];
    const rawLabels = Array.isArray(payload.detector_labels) ? payload.detector_labels.map(v => String(v)) : [];
    if (rawNames.length === 0 || rawLabels.length === 0) return detNames;

    const labelByName = new Map<string, string>();
    rawNames.forEach((name, idx) => {
        const label = rawLabels[idx];
        if (label && label !== 'NA') labelByName.set(name, label);
    });

    return detNames.map(det => labelByName.get(det) || det);
};

const App = ({ embedded = false, cockpitTheme = null }: { embedded?: boolean; cockpitTheme?: 'light' | 'dark' | null }) => {
    const [matrices, setMatrices] = useState<string[]>([]);
    const [currentFile, setCurrentFile] = useState('');
    const [sampleFiles, setSampleFiles] = useState<string[]>([]);
    const [currentSample, setCurrentSample] = useState('');
    const [matrix, setMatrix] = useState<MatrixRow[]>([]);
    const [loading, setLoading] = useState(true);
    const [detectors, setDetectors] = useState<string[]>([]);
    const [detectorLabels, setDetectorLabels] = useState<string[]>([]);
    const [unmixedData, setUnmixedData] = useState<DataRow[]>([]);
    const [rawData, setRawData] = useState<DataRow[]>([]);
    const [isUnmixingMatrix, setIsUnmixingMatrix] = useState(false);
    const [saveStatus, setSaveStatus] = useState<'idle' | 'saving' | 'saved'>('idle');
    const [errorMessage, setErrorMessage] = useState('');
    const [guiStateLoaded, setGuiStateLoaded] = useState(false);
    const [unmixingMethod, setUnmixingMethod] = useState<UnmixingMethod>('AutoSpectral');

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
    const unmixingMethodRef = useRef<UnmixingMethod>('AutoSpectral');

    useEffect(() => {
        const previousTitle = document.title;
        document.title = 'Spectreasy · Matrix Adjustment';
        return () => { document.title = previousTitle; };
    }, []);

    const fetchMatrices = async () => {
        const res = await axios.get(projectUrl('/matrices'));
        const list = Array.isArray(res.data) ? res.data : [];
        setMatrices(list);
        return list;
    };

    const fetchSamples = async () => {
        const res = await axios.get(projectUrl('/samples'));
        const list = Array.isArray(res.data) ? res.data : [];
        setSampleFiles(list);
        return list;
    };

    const fetchStatus = async () => {
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
        const resData = await axios.get(projectUrl(`/data${q}`));
        if (resData.data.error) {
            setRawData([]);
            setUnmixedData([]);
            return;
        }
        const raw = Array.isArray(resData.data.raw_data) ? resData.data.raw_data as DataRow[] : [];
        setRawData(raw);
        setDetectorLabels(alignDetectorLabels(detNames, resData.data || {}));
        if (resData.data.sample_name) setCurrentSample(asScalarString(resData.data.sample_name));
        await runUnmix(matrixData, raw, filenameForType);
    };

    const fetchData = async (filename = currentFile, sampleName = currentSample) => {
        if (!filename) {
            setLoading(false);
            return;
        }
        setLoading(true);
        setErrorMessage('');
        try {
            const resMatrix = await axios.get(projectUrl(`/load_matrix?filename=${encodeURIComponent(filename)}`));
            if (resMatrix.data?.error) throw new Error(asScalarString(resMatrix.data.error));
            const matrixData = Array.isArray(resMatrix.data) ? resMatrix.data as MatrixRow[] : [];
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
                await fetchStatus();
                const [mats, samples] = await Promise.all([fetchMatrices(), fetchSamples()]);
                await fetchUserGuiState();
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
        init();
        // Initial boot should run once; the called helpers intentionally read their current defaults.
        // eslint-disable-next-line react-hooks/exhaustive-deps
    }, []);

    async function runUnmix(currentM: MatrixRow[], currentRaw: DataRow[], filename: string | boolean = currentFile) {
        if (!Array.isArray(currentM) || currentM.length === 0 || !Array.isArray(currentRaw) || currentRaw.length === 0) {
            setUnmixedData([]);
            return false;
        }
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
            if (Array.isArray(res.data)) {
                setUnmixedData(res.data as DataRow[]);
                return true;
            }
            throw new Error(asScalarString(res.data?.error, 'The R backend returned no unmixed data.'));
        } catch (error) {
            setUnmixedData([]);
            setErrorMessage(error instanceof Error ? error.message : 'Sample unmixing failed.');
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
            await fetchMatrices();
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
    const chartWidth = Math.max(1040, detectors.length * 22);
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

    if (loading) return (
        <div style={{
            height: '100vh',
            width: '100%',
            maxWidth: '100%',
            display: 'flex',
            alignItems: 'center',
            justifyContent: 'center',
            backgroundColor: g.bgGradient,
            backgroundImage: renderedTheme === 'dark'
                ? 'linear-gradient(90deg, rgba(255,255,255,.02) 1px, transparent 1px), linear-gradient(rgba(255,255,255,.02) 1px, transparent 1px)'
                : 'linear-gradient(90deg, rgba(38,63,115,.035) 1px, transparent 1px), linear-gradient(rgba(38,63,115,.035) 1px, transparent 1px)',
            backgroundSize: '22px 22px',
            color: g.accent,
            fontFamily: 'Avenir Next, "Segoe UI", sans-serif',
            fontWeight: 500
        }}>
            <RefreshCw className="animate-spin" style={{ marginRight: 12 }} size={24} /> Initializing Workspace...
        </div>
    );

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
                            <div style={{ width: '100%', minWidth: 0 }}>
                                <svg
                                    ref={svgRef}
                                    style={{ display: 'block', width: '100%', maxWidth: '100%' }}
                                    width="100%"
                                    height={chartHeight + 56}
                                    viewBox={`0 0 ${chartWidth} ${chartHeight + 56}`}
                                >
                                    {[0, 25, 50, 75, 100].map(tick => {
                                        const y = chartHeight - (tick / 100) * (chartHeight - 32) - 24;
                                        return (
                                            <g key={tick}>
                                                <line x1={spectrumLeft} y1={y} x2={spectrumRight} y2={y} stroke={g.gridLine} strokeWidth={1} />
                                            </g>
                                        );
                                    })}
                                    {signatureDetectorLabels.map((label, idx) => {
                                        const x = detectorX(idx);
                                        return (
                                            <g key={`${label}-${idx}`}>
                                                <line x1={x} y1={6} x2={x} y2={chartHeight - 24} stroke={g.gridLine} strokeWidth={1} />
                                                <text x={x} y={chartHeight + 4} fontSize={11} textAnchor="end" transform={`rotate(-90 ${x} ${chartHeight + 4})`} fill={g.textMuted}>
                                                    {label}
                                                </text>
                                            </g>
                                        );
                                    })}
                                    <line x1={spectrumLeft} y1={chartHeight - 24} x2={spectrumRight} y2={chartHeight - 24} stroke={g.glassBorder} strokeWidth={3} />
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
                        <div style={{ ...glassCard, padding: 12, overflow: 'auto', maxHeight: 'none' }}>
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
                        </div>
                    </div>
                </div>
            </div>
        </div>
    );
};

export default App;
