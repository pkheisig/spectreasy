import type { AnalysisIdentitySignature } from './types'

const TEMPLATE_COLORS = ['#197783', '#d06d32', '#7d58a6', '#4f8f45', '#bc4d65', '#416fae']

type IdentityDefinition = {
  name: string
  positive: string[][]
  negative: string[][]
}
const COMMON_IMMUNE_IDENTITIES: IdentityDefinition[] = [
  {
    name: 'CD4 T cell',
    positive: [['CD3'], ['CD4']],
    negative: [['CD8'], ['CD19', 'CD20'], ['CD56'], ['CD14']],
  },
  {
    name: 'CD8 T cell',
    positive: [['CD3'], ['CD8']],
    negative: [['CD4'], ['CD19', 'CD20'], ['CD56'], ['CD14']],
  },
  {
    name: 'B cell',
    positive: [['CD19', 'CD20']],
    negative: [['CD3'], ['CD56'], ['CD14']],
  },
  {
    name: 'NK cell',
    positive: [['CD56']],
    negative: [['CD3'], ['CD19', 'CD20'], ['CD14']],
  },
  {
    name: 'Monocyte',
    positive: [['CD14', 'CD16']],
    negative: [['CD3'], ['CD19', 'CD20'], ['CD56']],
  },
]

function canonicalMarker(marker: string): string {
  return marker
    .toUpperCase()
    .replace(/\s+/g, '')
    .replace(/-(A|H|W)$/u, '')
    .replace(/[^A-Z0-9]/gu, '')
}

function firstAvailable(group: string[], available: Map<string, string>): string | undefined {
  for (const candidate of group) {
    const marker = available.get(canonicalMarker(candidate))
    if (marker) return marker
  }
  return undefined
}

export function commonImmuneIdentityTemplate(markers: string[]): AnalysisIdentitySignature[] {
  const available = new Map<string, string>()
  for (const marker of markers) {
    const canonical = canonicalMarker(marker)
    if (canonical && !available.has(canonical)) available.set(canonical, marker)
  }

  return COMMON_IMMUNE_IDENTITIES.flatMap((definition, index) => {
    const positive = definition.positive
      .map((group) => firstAvailable(group, available))
      .filter((marker): marker is string => Boolean(marker))
    if (positive.length !== definition.positive.length) return []
    const negative = definition.negative
      .map((group) => firstAvailable(group, available))
      .filter((marker): marker is string => Boolean(marker))
      .filter((marker) => !positive.includes(marker))
    return [{
      id: `immune-template-${canonicalMarker(definition.name).toLowerCase()}`,
      name: definition.name,
      color: TEMPLATE_COLORS[index % TEMPLATE_COLORS.length],
      positive_markers: positive,
      negative_markers: negative,
    }]
  })
}
