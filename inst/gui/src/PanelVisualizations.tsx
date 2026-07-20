import type { Dispatch, SetStateAction } from 'react';
import {
    bandColor,
    detectorColumnCenterX,
    detectorPointX,
    formatMetric,
    getSimilarityStyle,
    laserLabel,
    laserWavelength,
    linePath,
    signatureBandBins,
    signatureY,
    toNumber,
    toSimilarityValue,
    wavelengthToColor,
} from './panelBuilderShared';
import type { NumericRow, PanelPayload, TabId } from './panelBuilderShared';

export interface SelectedPanelEntry {
    fluor: string;
    slotIndex: number;
    marker: string;
    color: string;
    peakLaser: string;
    peakEmission: number;
}

interface PanelVisualizationsProps {
    payload: PanelPayload;
    selected: string[];
    selectedEntries: SelectedPanelEntry[];
    emissions: number[];
    lasers: string[];
    tab: TabId;
    setTab: Dispatch<SetStateAction<TabId>>;
    setMarkers: Dispatch<SetStateAction<Record<number, string>>>;
    spectraByName: Map<string, NumericRow>;
    similarityByName: Map<string, NumericRow>;
    colorByFluor: Map<string, string>;
    hoveredFluor: string | null;
    setHoveredFluor: Dispatch<SetStateAction<string | null>>;
    theme: 'light' | 'dark';
    error: string;
}

export function PanelVisualizations({
    payload,
    selected,
    selectedEntries,
    emissions,
    lasers,
    tab,
    setTab,
    setMarkers,
    spectraByName,
    similarityByName,
    colorByFluor,
    hoveredFluor,
    setHoveredFluor,
    theme,
    error,
}: PanelVisualizationsProps) {
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
    );
}
