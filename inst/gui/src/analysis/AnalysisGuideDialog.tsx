import { useEffect, useMemo, useRef, useState, type ReactNode } from 'react'
import { BookOpen, ChevronRight, X } from 'lucide-react'
import { createPortal } from 'react-dom'
import type { AnalysisMethod } from './types'

type GuideSection = {
  id: string
  label: string
  title: string
  content: ReactNode
}

export function AnalysisGuideDialog({
  methods,
  onClose,
}: {
  methods: AnalysisMethod[]
  onClose: () => void
}) {
  const [activeId, setActiveId] = useState('workflow')
  const dialogRef = useRef<HTMLElement>(null)
  const onCloseRef = useRef(onClose)
  const sections = useMemo<GuideSection[]>(() => {
    const visibleMethods = methods.filter((method) => method.visible !== false)
    const ready = visibleMethods.filter((method) => method.available)
    const unavailable = visibleMethods.filter((method) => !method.available)
    return [
      {
        id: 'gating',
        label: 'Gating workspace',
        title: 'Build the population before analyzing it',
        content: <>
          <p>Open this workspace from <strong>Other tools → Population analysis</strong>. Gating, downstream analysis, cell identities, statistics, and exports remain together here.</p>
          <p>The left sidebar is the population hierarchy. Selecting a population filters every plot assigned to it and makes new gates children of that population.</p>
          <ol>
            <li>Add as many fixed-size square plots as needed.</li>
            <li>For each plot, choose its population, axes, linear, asinh, or biexponential transforms, and display type.</li>
            <li>Use a density or continuously marker-colored scatter, density contours, hexbin, or a one-marker histogram.</li>
            <li>Draw rectangle, ellipse, or polygon gates on a 2D plot, or a range gate on a histogram, then name the child population.</li>
          </ol>
          <p>Click the underlined X or Y dimension outside the plot canvas to choose another channel, matching the control-gating interaction. Scatter plots use a clean, gridless canvas; numeric ticks remain available around the data region.</p>
          <p>Type in the dimension menu to filter long panels by either detector channel or marker name. Undo and Redo restore complete workspace changes; use Ctrl/Cmd+Z and Ctrl/Cmd+Shift+Z outside text fields.</p>
          <p><strong>Save</strong> and <strong>Load</strong> beside Populations handle reusable gate CSVs. They contain only the nested population hierarchy, semantic roles, channels, optional file scope, and geometry. Loading checks the selected source first, shows the exact replacement and compatibility adjustments, and replaces gates only after confirmation. Undo restores the previous tree during the session.</p>
          <p>Delete any gating or result plot with its trash button. Deleting every plot does not delete the gates or analysis object; use Add plot to create another view. Filter or collapse the population hierarchy when it grows. Exact gate coordinates can be edited in the Population inspector; histogram range boundaries can also be dragged directly.</p>
          <p><strong>Backgating</strong> overlays another population without changing the events used by the active plot. Gates can apply to every compatible file or only the current file. Plots can be duplicated and reordered without changing their data.</p>
          <p className="analysis-guide-callout">The highlighted plot is the one edited by Plot settings. Its border replaces a separate “active plot” title.</p>
        </>,
      },
      {
        id: 'workflow',
        label: 'Analysis workflow',
        title: 'Gate first, analyze a population second',
        content: <>
          <p>The analysis workspace keeps manual gating separate from dimensional reduction, clustering, and trajectory inference.</p>
          <ol>
            <li>Create or select a population in the hierarchy.</li>
            <li>Choose <strong>Analyze population</strong>.</li>
            <li>Select one or more compatible FCS files to pool, select markers, then choose a clustering method from the dropdown or check <strong>Skip</strong>.</li>
            <li>Choose the dimensional-reduction map independently and run the pipeline.</li>
          </ol>
          <p>Open <strong>Advanced settings</strong> to tune only the selected algorithms. Each field is constrained to values accepted by its backend; Reset restores that method’s maintained defaults.</p>
          <p>Pooled runs apply the same population gate to each file before seeded sampling. Every result row retains its source file and original FCS event row.</p>
          <p className="analysis-guide-callout">Clustering and embedding objects are saved separately. Keep the clustering selection and choose another map to reuse the same clustering.</p>
        </>,
      },
      {
        id: 'settings',
        label: 'Advanced settings',
        title: 'Tune only the methods you selected',
        content: <>
          <p>Open <strong>Advanced settings</strong> below the pipeline selectors. The panel changes with the chosen clustering, map, or trajectory method, so irrelevant parameters are never shown.</p>
          <p>Every setting has a maintained default and an allowed type, range, or choice. Invalid combinations—such as a final FlowSOM learning rate above its initial rate—disable Run and explain what must be corrected.</p>
          <p><strong>Reset</strong> restores only that method. Changing a setting produces a distinct cached object; returning to a previously completed configuration reuses its saved object.</p>
          <p>The common Run settings control the maximum analyzed events, arcsinh cofactor, and master seed. Method-specific values such as neighbors, clusters, iterations, distance metrics, or diffusion components live only in Advanced settings.</p>
        </>,
      },
      {
        id: 'availability',
        label: 'Method availability',
        title: 'Only executable methods can be selected',
        content: <>
          <p><strong>Ready</strong> means the adapter, required runtime, and small-population execution contract have all passed. These methods can be selected and run.</p>
          <p><strong>Setup required</strong> means the adapter is supported but a runtime dependency must be installed. <strong>Unavailable</strong> means the method is not implemented in this build. Both states are gray and cannot be clicked.</p>
          <h4>Ready in this runtime</h4>
          <ul>{ready.map((method) => <li key={method.id}><strong>{method.name}</strong></li>)}</ul>
          {unavailable.length ? <><h4>Currently unavailable</h4><ul>{unavailable.map((method) => <li key={method.id}><strong>{method.name}</strong> — {method.blocker || method.next_action}</li>)}</ul></> : null}
        </>,
      },
      {
        id: 'prerequisites',
        label: 'Prerequisites',
        title: 'Automatic stages and user actions are different',
        content: <>
          <p>Clustering and dimensional reduction both use the same sampled, transformed marker matrix. Cluster labels are joined to the selected map by event ID; the embedding is not calculated from the cluster labels.</p>
          <p>Pipeline prerequisites such as <strong>diffusion map → DPT</strong> are automatic. Spectreasy creates the intermediate model, sends that exact object into the next stage, and records the hand-off. Methods that also require clustering show explicit cluster and map selectors.</p>
          <p>For trajectory analysis, choose <strong>Set root</strong> in the gating toolbar and click an event inside the population being analyzed. The Run button remains disabled until the root and all required stages are valid.</p>
          <p>Slingshot and TSCAN ask for both clustering and dimensional-reduction inputs. DPT automatically creates and retains its diffusion map. Other trajectory methods build or reuse the neighbor graph, diffusion representation, or cluster input required by their backend.</p>
          <p>Slingshot uses robust Euclidean cluster-center distances by default; covariance-aware distance remains available in Advanced settings. If a valid PAGA input has no inter-cluster kNN edges, Spectreasy records a centroid minimum-spanning-tree fallback and still runs DPT on the original event-level diffusion graph.</p>
        </>,
      },
      {
        id: 'visualization',
        label: 'Result plots',
        title: 'Coordinates and colors are independent',
        content: <>
          <p>Result plots start in the compact square 2D view. Open the coordinate control to assign X and Y. When an analysis returns at least three coordinates, switch to 3D there and assign Z.</p>
          <p>The 3D view uses Plotly for rotation, zooming, panning, and hover inspection. Methods or marker selections that produce only two coordinates never show a 3D control.</p>
          <p>Color can show density, clusters, pseudotime, a measured marker, or predicted cell identity when those values exist. Continuous colors have an independent palette picker.</p>
          <p>Available continuous palettes include the control-gating gradient, Viridis, Sunset, Plasma, Inferno, Magma, Cividis, Turbo, Ice–fire, and Spectral.</p>
          <p>Use <strong>Add plot</strong> for another independent view of the same result. Each plot can use different coordinates, dimensionality, color, and palette without rerunning the analysis.</p>
          <p>Export a plot as PNG, vector SVG, or interactive HTML. <strong>Export plot data</strong> writes the complete event table with coordinates, marker values, cluster and trajectory assignments, and cell identities when present.</p>
        </>,
      },
      {
        id: 'identities',
        label: 'Cell identities',
        title: 'Annotate with panel-specific marker evidence',
        content: <>
          <p>The Cell identities tab becomes available only after a map or trajectory result exists. Labels are calculated from the transformed marker values—not from the visual map coordinates.</p>
          <ol>
            <li>If clustering was run, use <strong>Discover cluster markers</strong> to rank cluster-versus-rest marker evidence before naming anything.</li>
            <li>Load the panel-aware common immune template or name at least two candidate identities manually. Templates are editable and identities without the required positive markers are omitted.</li>
            <li>Choose markers expected to be high under <strong>POS</strong> and exclusion markers expected to be low under <strong>NEG</strong>.</li>
            <li>Open the gear only when the default assignment strictness needs adjustment.</li>
            <li>Run annotation, then color the result by <strong>Predicted identity</strong> or by individual markers to verify it.</li>
          </ol>
          <div className="analysis-guide-definitions">
            <div><strong>Minimum marker separation score</strong><p>Cluster-versus-rest rank AUC. 0.5 means no separation; 1 means complete separation. Marker discovery explains clusters but does not itself name cell types.</p></div>
            <div><strong>Minimum identity match</strong><p>How well the cell must fit its best marker pattern. Higher values make assignment more selective.</p></div>
            <div><strong>Minimum lead over second choice</strong><p>How far the best identity must beat the runner-up. Higher values protect against ambiguous labels.</p></div>
            <div><strong>Marker sensitivity</strong><p>How strongly expression above or below the population’s robust center contributes. 1 is balanced; higher values make smaller marker differences count more strongly.</p></div>
          </div>
          <p className="analysis-guide-callout">A cell must pass both confidence checks. Cells that do not remain Unassigned; the map coordinates themselves never influence identity.</p>
        </>,
      },
      {
        id: 'statistics',
        label: 'Statistics & export',
        title: 'Summarize, compare, and export populations',
        content: <>
          <p>The Population inspector reports event count, percentage of the parent and total population, plus marker median, mean, and robust standard deviation.</p>
          <p>For a staining index, create sibling gates under the same parent and assign one the <strong>POS</strong> role and one the <strong>NEG</strong> role. Spectreasy calculates:</p>
          <p className="analysis-guide-formula">(median POS − median NEG) / (2 × robust SD of NEG)</p>
          <p>The Export tab writes the active population as FCS, CSV, or both. Set Maximum events to 0 for every event or use seeded random sampling without replacement. The same population can be exported from one file or every file in the selected source. Use the folder button to choose a destination inside the active project.</p>
          <p>Hierarchy-statistics CSVs include every population and marker summary. Event exports record the source, original event row, population path, seed, software writer, and checksum.</p>
          <p><strong>Gate CSV</strong> is the ordinary sharing format for reusable gating templates. It deliberately excludes plots, annotations, active selections, analysis settings, seeds, and trajectory-root state. Empty source-file cells are global gates; populated values contain only an FCS filename. Spectreasy resolves that name inside the selected source or requires a unique match in the active project.</p>
          <p><strong>Export workspace JSON</strong> is the advanced whole-workspace option. It preserves the complete hierarchy, gates, roles, plots, transforms, annotations, active selections, settings, seeds, and trajectory root. The normal internal state continues to autosave at <strong>.spectreasy/analysis-v2/workspace.json</strong>.</p>
        </>,
      },
      {
        id: 'reproducibility',
        label: 'Reproducibility',
        title: 'Seeds and intermediate outputs are preserved',
        content: <>
          <p>Downsampling is random without replacement and is derived from the master seed, file, population, markers, and transform. Changing only the map therefore keeps the same sampled events.</p>
          <p>Stochastic algorithms receive the recorded seed. Fitted clustering, embedding, and trajectory objects, event-level values, advanced settings, and DPT diffusion coordinates are saved with their parameters and checksums.</p>
        </>,
      },
      {
        id: 'code',
        label: 'R code access',
        title: 'The GUI and R use the same backend',
        content: <>
          <p>The browser does not contain a separate scientific implementation. Method discovery, clustering, dimensional reduction, trajectory inference, annotation, validation, caching, and provenance call the same Spectreasy package backend available from R.</p>
          <ul>
            <li><strong>analysis_workspace()</strong>, <strong>add_population_gate()</strong>, <strong>update_population_gate()</strong>, and <strong>delete_population_gate()</strong> provide gate-hierarchy CRUD.</li>
            <li><strong>export_population_gates()</strong> and <strong>import_population_gates()</strong> save and replace the same schema-versioned gate CSV used by the sidebar controls.</li>
            <li><strong>population_statistics()</strong>, <strong>staining_index()</strong>, and <strong>export_gated_population()</strong> reproduce the gating inspector and export actions.</li>
            <li><strong>analysis_methods()</strong> reports availability and every advanced-parameter schema.</li>
            <li><strong>install_analysis_dependencies()</strong> installs all optional R methods and the pinned private Python runtime; <strong>analysis_runtime_status()</strong> inspects that runtime without changing system Python.</li>
            <li><strong>analyze_population()</strong> runs clustering, maps, and trajectories with named setting lists.</li>
            <li><strong>find_population_markers()</strong> ranks cluster-versus-rest marker evidence and exports the complete table.</li>
            <li><strong>population_identity_templates()</strong> creates the same panel-aware, editable common immune starting patterns used by the GUI.</li>
            <li><strong>annotate_population()</strong> applies the same marker-pattern identity model.</li>
            <li><strong>plot_population_analysis()</strong> and <strong>export_population_analysis()</strong> reproduce plot customization and table or figure export in code.</li>
          </ul>
          <p>Saved results can be reopened with <strong>load_population_analysis()</strong>. Interactive plot rendering remains in the browser, while all scientific result values originate from the shared backend.</p>
        </>,
      },
    ]
  }, [methods])
  const active = sections.find((section) => section.id === activeId) ?? sections[0]

  useEffect(() => {
    onCloseRef.current = onClose
  }, [onClose])

  useEffect(() => {
    dialogRef.current?.querySelector<HTMLButtonElement>('[role="tab"]')?.focus()
    const keydown = (event: KeyboardEvent) => {
      if (event.key === 'Escape') {
        event.preventDefault()
        onCloseRef.current()
        return
      }
      if (event.key !== 'Tab') return
      const focusable = Array.from(dialogRef.current?.querySelectorAll<HTMLElement>('button:not([disabled]), [href], [tabindex]:not([tabindex="-1"])') ?? [])
      if (!focusable.length) return
      const first = focusable[0]
      const last = focusable[focusable.length - 1]
      if (event.shiftKey && document.activeElement === first) {
        event.preventDefault()
        last.focus()
      } else if (!event.shiftKey && document.activeElement === last) {
        event.preventDefault()
        first.focus()
      }
    }
    document.addEventListener('keydown', keydown)
    return () => document.removeEventListener('keydown', keydown)
  }, [])

  return createPortal(
    <div className="analysis-guide-overlay" role="presentation" onMouseDown={onClose}>
      <section ref={dialogRef} className="analysis-guide-dialog" role="dialog" aria-modal="true" aria-labelledby="analysis-guide-title" onMouseDown={(event) => event.stopPropagation()}>
        <header>
          <div><BookOpen size={20} /><span><small>Population analysis</small><strong id="analysis-guide-title">Guide</strong></span></div>
          <button type="button" onClick={onClose} aria-label="Close analysis guide"><X size={18} /></button>
        </header>
        <div className="analysis-guide-layout">
          <nav role="tablist" aria-label="Analysis guide sections">
            {sections.map((section) => (
              <button
                key={section.id}
                type="button"
                role="tab"
                aria-selected={section.id === active.id}
                tabIndex={section.id === active.id ? 0 : -1}
                className={section.id === active.id ? 'is-active' : ''}
                onClick={(event) => {
                  setActiveId(section.id)
                  event.currentTarget.focus()
                }}
              >
                <span>{section.label}</span><ChevronRight size={14} />
              </button>
            ))}
          </nav>
          <article role="tabpanel" tabIndex={0}>
            <small>{String(sections.findIndex((section) => section.id === active.id) + 1).padStart(2, '0')}</small>
            <h3>{active.title}</h3>
            {active.content}
          </article>
        </div>
      </section>
    </div>,
    document.body,
  )
}
