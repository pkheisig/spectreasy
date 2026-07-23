import {
  compactDetectorLabel,
  DETECTOR_LABEL_ROTATION,
  detectorLaserSegments,
  detectorSpectralColors,
} from './detectorAxis'
import type { DetectorAxisEntry } from './detectorAxis'

type DetectorSpectrumAxisProps = {
  entries: DetectorAxisEntry[]
  cytometer?: unknown
  xForIndex: (index: number) => number
  left: number
  right: number
  plotTop: number
  baselineY: number
  gridColor: string
  axisColor: string
  textColor: string
  gradientId: string
}

export function DetectorSpectrumAxis({
  entries,
  cytometer = '',
  xForIndex,
  left,
  right,
  plotTop,
  baselineY,
  gridColor,
  axisColor,
  textColor,
  gradientId,
}: DetectorSpectrumAxisProps) {
  if (entries.length === 0) return null
  const detectorColors = detectorSpectralColors(entries)
  const laserSegments = detectorLaserSegments(entries, cytometer)
  const detectorStep = entries.length <= 1 ? right - left : xForIndex(1) - xForIndex(0)
  const spectrumY = baselineY + 8
  const laserY = spectrumY + 10
  const laserHeight = 18
  const labelY = laserY + laserHeight + 8
  const segmentStartX = (index: number) => index === 0
    ? left
    : Math.max(left, xForIndex(index) - detectorStep / 2)
  const segmentEndX = (index: number) => index === entries.length - 1
    ? right
    : Math.min(right, xForIndex(index) + detectorStep / 2)

  return <g className="detector-spectrum-axis">
    <defs>
      <linearGradient id={gradientId} x1="0%" y1="0%" x2="100%" y2="0%">
        {detectorColors.map((color, index) => {
          const percentage = (index / Math.max(1, entries.length - 1)) * 100
          return <stop key={`${entries[index].detector}-${index}`} offset={`${percentage.toFixed(2)}%`} stopColor={color} />
        })}
      </linearGradient>
    </defs>
    {entries.map((entry, index) => {
      const x = xForIndex(index)
      return <g key={`${entry.detector}-${index}`}>
        <line x1={x} y1={plotTop} x2={x} y2={baselineY} stroke={gridColor} strokeWidth={1} />
        <text
          x={x}
          y={labelY}
          fontSize={9.5}
          fontWeight={650}
          textAnchor="end"
          dominantBaseline="middle"
          transform={`rotate(${DETECTOR_LABEL_ROTATION} ${x} ${labelY})`}
          fill={textColor}
        >
          {compactDetectorLabel(entry.label || entry.detector)}
        </text>
      </g>
    })}
    <line x1={left} y1={baselineY} x2={right} y2={baselineY} stroke={axisColor} strokeWidth={3} />
    <rect x={left} y={spectrumY} width={right - left} height={7} rx={1.5} fill={`url(#${gradientId})`} />
    {laserSegments.map((segment) => {
      const x1 = segmentStartX(segment.startIndex)
      const x2 = segmentEndX(segment.endIndex)
      const width = Math.max(1, x2 - x1)
      return <g key={`${segment.key}-${segment.startIndex}`}>
        {segment.startIndex > 0 && <line x1={x1} y1={plotTop} x2={x1} y2={laserY + laserHeight} stroke={segment.color} strokeWidth={1.5} strokeOpacity={0.52} />}
        <rect x={x1 + 1} y={laserY} width={Math.max(1, width - 2)} height={laserHeight} rx={3} fill={segment.color} />
        <text x={x1 + width / 2} y={laserY + 12.5} textAnchor="middle" fontSize={10.5} fontWeight={800} fill={segment.textColor}>
          {segment.label} <tspan fontWeight={600} opacity={0.82}>{segment.wavelength}</tspan>
        </text>
      </g>
    })}
  </g>
}
