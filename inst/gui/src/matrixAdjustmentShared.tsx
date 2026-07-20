/* eslint-disable react-refresh/only-export-components -- shared applet primitives intentionally colocate one dropdown with its data helpers */
import { useState, useEffect, useId, useRef } from 'react';
import { Check, ChevronDown } from 'lucide-react';
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

export {
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
};
export type { DataRow, MatrixPayload, MatrixRow, UnmixingMethod, ViewConfig };
