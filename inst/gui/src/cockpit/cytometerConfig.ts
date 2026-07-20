export const cytometerOptions = [
  ["auto", "Auto"],
  ["aurora", "Cytek Aurora"],
  ["northern_lights", "Cytek Northern Lights"],
  ["id7000", "Sony ID7000"],
  ["discover_s8", "BD FACSDiscover S8"],
  ["discover_a8", "BD FACSDiscover A8"],
  ["a5se", "BD FACSymphony A5 SE"],
  ["opteon", "Agilent NovoCyte Opteon"],
  ["mosaic", "Beckman Coulter CytoFLEX Mosaic"],
  ["xenith", "Thermo Fisher Attune Xenith"],
] as const;

export const cytometerLabels = Object.fromEntries(cytometerOptions) as Record<string, string>;

export function normalizeCockpitCytometer(value: unknown): string {
  const token = String(value ?? "auto")
    .trim()
    .toLowerCase()
    .replace(/[^a-z0-9]+/g, "_")
    .replace(/^_|_$/g, "");
  const legacyAliases: Record<string, string> = {
    aurora_5l: "aurora",
    aurora_4l: "aurora",
    cytek_aurora: "aurora",
    discover: "auto",
    cytek_aurora_discover: "auto",
    bd_facsdiscover_xenith: "xenith",
    thermo_fisher_attune_xenith: "xenith",
  };
  const normalized = legacyAliases[token] ?? token;
  return normalized in cytometerLabels ? normalized : "auto";
}
