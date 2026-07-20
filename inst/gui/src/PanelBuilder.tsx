import { useEffect, useMemo, useRef, useState } from 'react';
import type { CSSProperties, PointerEvent as ReactPointerEvent } from 'react';
import axios from 'axios';
import { Moon, PanelLeftClose, PanelLeftOpen, Plus, Save, Sun, Trash2, Upload } from 'lucide-react';
import './PanelBuilder.css';
import { PanelVisualizations } from './PanelVisualizations';
import {
    API_BASE,
    PdfIcon,
    base64ToBlob,
    binEmission,
    csvEscape,
    detectImportedPanelRows,
    downloadBlob,
    emptySlots,
    getCytometerName,
    laserOrder,
    mapDetectorToEmission,
    normalizeMarkers,
    panelProjectBody,
    panelProjectUrl,
    unique,
    unboxGuiState,
} from './panelBuilderShared';
import type {
    FluorInfo,
    NumericRow,
    PanelExportResponse,
    PanelPayload,
    TabId,
} from './panelBuilderShared';

const PanelBuilder = ({ embedded = false, cockpitTheme = null }: { embedded?: boolean; cockpitTheme?: 'light' | 'dark' | null }) => {
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
        if (embedded && cockpitTheme) return cockpitTheme;
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
        if (embedded) return;
        localStorage.setItem('spectreasy-theme', theme);
        localStorage.removeItem('spectreasy_theme');
        document.documentElement.dataset.theme = theme;
    }, [embedded, theme]);

    useEffect(() => {
        localStorage.setItem('spectreasy_slots', JSON.stringify(slots));
    }, [slots]);

    useEffect(() => {
        localStorage.setItem('spectreasy_markers', JSON.stringify(markers));
    }, [markers]);

    useEffect(() => {
        if (!guiStateLoaded) return;
        const timer = window.setTimeout(() => {
            void axios.post(`${API_BASE}/gui_state`, panelProjectBody({
                module: 'panel_builder',
                config_json: {
                    cytometer: getCytometerName(cytometer),
                    configuration: getCytometerName(configuration),
                    ...(!embedded ? { theme } : {}),
                    slots,
                    markers,
                    tab,
                    sidebarWidth,
                    sidebarCollapsed
                }
            })).catch(() => null);
        }, 500);
        return () => window.clearTimeout(timer);
    }, [cytometer, configuration, theme, slots, markers, tab, sidebarWidth, sidebarCollapsed, guiStateLoaded, embedded]);

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
        const res = await axios.post(`${API_BASE}/spectral_panel_metrics`, panelProjectBody({
            cytometer: nextCytometer,
            configuration: nextConfiguration,
            fluorophores: nextSelected,
        }));
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
                const stateRes = await axios.get(panelProjectUrl('/gui_state?module=panel_builder')).catch(() => null);
                const saved = unboxGuiState(stateRes?.data?.config || {}) as Record<string, unknown>;
                const defaults = bootDefaultsRef.current;
                const savedCytometer = typeof saved.cytometer === 'string' ? getCytometerName(saved.cytometer) : defaults.cytometer;
                const savedConfiguration = typeof saved.configuration === 'string' ? getCytometerName(saved.configuration) : defaults.configuration;
                const savedSlots = Array.isArray(saved.slots) ? saved.slots.map(String) : defaults.slots;
                const savedMarkers = saved.markers && typeof saved.markers === 'object' ? normalizeMarkers(saved.markers) : defaults.markers;
                if (!embedded && (saved.theme === 'light' || saved.theme === 'dark')) setTheme(saved.theme);
                if (saved.tab === 'panel' || saved.tab === 'similarity' || saved.tab === 'signatures') setTab(saved.tab);
                if (typeof saved.sidebarWidth === 'number' && Number.isFinite(saved.sidebarWidth)) {
                    setSidebarWidth(Math.min(440, Math.max(180, saved.sidebarWidth)));
                }
                if (typeof saved.sidebarCollapsed === 'boolean') setSidebarCollapsed(saved.sidebarCollapsed);
                setSlots(savedSlots);
                slotsRef.current = savedSlots;
                setMarkers(savedMarkers);
                const res = await axios.post(`${API_BASE}/spectral_panel_metrics`, panelProjectBody({
                    cytometer: savedCytometer,
                    configuration: savedConfiguration,
                    fluorophores: savedSlots.filter(Boolean),
                }));
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
    }, [embedded]);

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
            const res = await axios.post(`${API_BASE}/export_spectral_panel_overview`, panelProjectBody({
                cytometer,
                configuration,
                fluorophores: selectedRows.map(row => row.fluor),
                markers: selectedRows.map(row => row.marker),
            }));
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

    return (
        <div className={`panel-builder ${embedded && cockpitTheme ? cockpitTheme : theme}`}>
            <header className="panel-topbar">
                <div>
                    <h1>Spectral Panel Builder</h1>
                    <p>{selected.length} fluorophores selected{selectedConfigurationLabel ? ` / ${selectedConfigurationLabel}` : ''}</p>
                </div>
                <div className="panel-actions">
                    {!embedded && <button
                        type="button"
                        className="export-button"
                        onClick={() => setTheme(prev => prev === 'light' ? 'dark' : 'light')}
                        aria-label="Toggle theme"
                        style={{ padding: '0 10px', width: '40px' }}
                    >
                        {theme === 'light' ? <Moon size={16} /> : <Sun size={16} />}
                    </button>}
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

                <PanelVisualizations
                    payload={payload}
                    selected={selected}
                    selectedEntries={selectedEntries}
                    emissions={emissions}
                    lasers={lasers}
                    tab={tab}
                    setTab={setTab}
                    setMarkers={setMarkers}
                    spectraByName={spectraByName}
                    similarityByName={similarityByName}
                    colorByFluor={colorByFluor}
                    hoveredFluor={hoveredFluor}
                    setHoveredFluor={setHoveredFluor}
                    theme={embedded && cockpitTheme ? cockpitTheme : theme}
                    error={error}
                />
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
