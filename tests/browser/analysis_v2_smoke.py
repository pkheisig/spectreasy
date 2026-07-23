import sys
from pathlib import Path

from playwright.sync_api import sync_playwright


port = int(sys.argv[1])
captures = [Path(value) for value in sys.argv[2:]]
errors: list[str] = []


def capture(page, index: int) -> None:
    if index < len(captures):
        captures[index].parent.mkdir(parents=True, exist_ok=True)
        page.screenshot(path=str(captures[index]), full_page=True)


with sync_playwright() as playwright:
    browser = playwright.chromium.launch(headless=True)
    page = browser.new_page(viewport={"width": 1600, "height": 1000}, device_scale_factor=1)
    page.on("console", lambda message: errors.append(message.text) if message.type == "error" else None)
    page.on("pageerror", lambda error: errors.append(str(error)))
    page.goto(
        f"http://127.0.0.1:{port}/?api=http%3A%2F%2F127.0.0.1%3A{port}#token=analysis-v2-smoke-token",
        wait_until="networkidle",
    )

    project_prompt = page.locator(".project-initialize-confirm")
    if project_prompt.count() and project_prompt.is_visible():
        project_prompt.locator(".button-ghost", has_text="Not now").click()
    page.get_by_role("button", name="Samples").click()
    assert page.get_by_role("button", name="Analyze samples").count() == 0
    page.get_by_role("button", name="Other tools").click()
    population_entry = page.locator(".rail-subitem", has_text="Population analysis")
    assert population_entry.is_visible()
    assert page.locator(".rail-subitem", has_text="Panel builder").is_visible()
    population_entry.click()
    page.get_by_text("v2 workspace", exact=True).wait_for()

    analyze_button = page.get_by_role("button", name="Analyze population", exact=True)
    analyze_button.wait_for(timeout=60000)
    plot_card = page.locator(".analysis-plot-card").first
    x_axis = plot_card.get_by_role("button", name="X axis", exact=False)
    y_axis = plot_card.get_by_role("button", name="Y axis", exact=False)
    assert "FSC" in x_axis.get_attribute("aria-label")
    assert "SSC" in y_axis.get_attribute("aria-label")
    x_axis.click()
    axis_menu = plot_card.get_by_role("listbox", name="Choose X axis")
    axis_menu.wait_for()
    capture(page, 9)
    axis_menu.get_by_role("option", name="CD3 · CD3-A", exact=True).click()
    assert "CD3" in x_axis.get_attribute("aria-label")
    plot_card.get_by_role("button", name="Delete plot", exact=True).click()
    page.get_by_text("No plots", exact=True).wait_for()
    assert page.locator(".analysis-plot-card").count() == 0
    page.get_by_role("button", name="Add plot", exact=True).first.click()
    assert page.locator(".analysis-plot-card").count() == 1
    analyze_button.click()
    dialog = page.locator(".analysis-method-dialog")
    dialog.wait_for()
    dialog.get_by_role("button", name="Guide", exact=True).click()
    guide = page.locator(".analysis-guide-dialog")
    guide.wait_for()
    for section in (
        "Gating workspace", "Analysis workflow", "Advanced settings",
        "Method availability", "Prerequisites", "Result plots",
        "Cell identities", "Statistics & export", "Reproducibility", "R code access",
    ):
        assert guide.get_by_role("tab", name=section, exact=True).is_visible(), section
    identity_guide_tab = guide.get_by_role("tab", name="Cell identities", exact=True)
    identity_guide_tab.click()
    assert identity_guide_tab.get_attribute("aria-selected") == "true"
    assert "is-active" in (identity_guide_tab.get_attribute("class") or "")
    workflow_guide_tab = guide.get_by_role("tab", name="Analysis workflow", exact=True)
    assert workflow_guide_tab.get_attribute("aria-selected") == "false"
    assert "is-active" not in (workflow_guide_tab.get_attribute("class") or "")
    assert identity_guide_tab.evaluate("(element) => document.activeElement === element")
    guide.get_by_role("tabpanel").hover()
    page.wait_for_timeout(100)
    assert guide.get_by_text("Minimum identity match", exact=True).is_visible()
    assert guide.get_by_text("Minimum lead over second choice", exact=True).is_visible()
    assert guide.get_by_text("Marker sensitivity", exact=True).is_visible()
    capture(page, 0)
    guide.get_by_role("button", name="Close analysis guide").click()
    file_picker = dialog.locator(".analysis-analysis-file-picker")
    assert file_picker.locator("label").count() == 2
    unchecked_files = [checkbox for checkbox in file_picker.get_by_role("checkbox").all() if not checkbox.is_checked()]
    assert len(unchecked_files) == 1
    unchecked_files[0].check()
    assert file_picker.get_by_role("checkbox").count() == 2
    assert all(checkbox.is_checked() for checkbox in file_picker.get_by_role("checkbox").all())
    cluster_stage = dialog.locator(".analysis-pipeline-stage").nth(0)
    skip = cluster_stage.get_by_role("checkbox", name="Skip")
    assert skip.is_visible()
    assert cluster_stage.get_by_role("option", name="Skip clustering").count() == 0

    cluster_select = dialog.get_by_label("Clustering method")
    map_select = dialog.get_by_label("Dimensional-reduction method")
    cluster_select.select_option("flowsom")
    map_select.select_option("umap")

    # Maintained dimensional-reduction adapters are selectable, including true HSNE.
    for method_id in ("pca", "tsne", "umap", "diffusion-map", "phate", "hsne"):
        assert map_select.locator(f'option[value="{method_id}"]').is_enabled(), method_id
    assert map_select.locator('option[value="tsne"]').inner_text() == "t-SNE"
    assert map_select.locator('option[value="hsne"]').inner_text() == "HSNE"
    assert not any("Pearson" in text or "McInnes" in text for text in map_select.locator("option").all_inner_texts())
    advanced = dialog.locator(".analysis-advanced-settings")
    advanced.locator("summary").click()
    assert advanced.locator(".analysis-advanced-method", has_text="FlowSOM").get_by_text("SOM grid width", exact=True).is_visible()
    umap_settings = advanced.locator(".analysis-advanced-method", has_text="UMAP")
    assert umap_settings.get_by_text("Minimum distance", exact=True).is_visible()
    umap_settings.get_by_text("Minimum distance", exact=True).locator("..").get_by_role("spinbutton").fill("0.08")
    capture(page, 1)
    advanced.locator("summary").click()

    dialog.get_by_role("button", name="Run pipeline").click()
    dialog.get_by_text("FlowSOM → UMAP", exact=True).wait_for(timeout=120000)
    result_plot = dialog.locator(".analysis-result-plot")
    result_plot.wait_for()
    result_box = result_plot.bounding_box()
    assert result_box is not None
    assert abs(result_box["width"] - result_box["height"]) < 2
    assert 360 <= result_box["width"] <= 380
    assert dialog.locator(".analysis-result-tabs").get_by_role("button", name="Cell identities").is_enabled()
    assert "spectreasy_outputs/analysis/" not in dialog.inner_text()

    color_select = dialog.get_by_label("Result color")
    color_labels = color_select.locator("option").all_inner_texts()
    assert "CD3" in color_labels
    assert "CD3-A" not in color_labels
    color_select.select_option(label="CD3")
    palette_select = dialog.get_by_label("Result color palette")
    assert {"Viridis", "Sunset"}.issubset(set(palette_select.locator("option").all_inner_texts()))
    palette_select.select_option("sunset")
    capture(page, 2)

    dialog.get_by_role("button", name="Add plot", exact=True).click()
    plot_cards = dialog.locator(".analysis-result-plot-card")
    assert plot_cards.count() == 2
    plot_cards.nth(1).get_by_role("button", name="Delete result plot", exact=True).click()
    assert plot_cards.count() == 1
    dialog.get_by_role("button", name="Add plot", exact=True).click()
    assert plot_cards.count() == 2
    for card in plot_cards.all():
        box = card.locator(".analysis-result-plot").bounding_box()
        assert box is not None and abs(box["width"] - box["height"]) < 2 and 360 <= box["width"] <= 380
    with page.expect_download() as table_download:
        dialog.get_by_role("button", name="Export plot data", exact=True).click()
    assert table_download.value.suggested_filename.endswith("-plot-data.csv")
    exported_table = Path(table_download.value.path()).read_text()
    assert "source_file" in exported_table
    assert "PKH_browser_sample.fcs" in exported_table
    assert "PKH_browser_sample_repeat.fcs" in exported_table
    with page.expect_download() as svg_download:
        plot_cards.nth(0).get_by_role("button", name="Export", exact=True).click()
        plot_cards.nth(0).get_by_role("button", name="SVG vector", exact=True).click()
    assert svg_download.value.suggested_filename.endswith(".svg")
    capture(page, 3)

    coordinate_picker = plot_cards.nth(0).get_by_label("Choose result coordinates")
    coordinate_picker.click()
    plot_cards.nth(0).get_by_label("Result Y coordinate").select_option("dimension_3")
    plot_cards.nth(0).get_by_label("Result X coordinate").select_option("dimension_2")
    plot_cards.nth(0).get_by_label("Result Y coordinate").select_option("dimension_1")
    assert plot_cards.nth(0).get_by_role("button", name="3D", exact=True).is_visible()
    plot_cards.nth(0).get_by_role("button", name="3D", exact=True).click()
    coordinate_picker.click()
    plot_cards.nth(0).locator(".analysis-result-plotly").wait_for(timeout=60000)
    plot_cards.nth(0).locator(".analysis-result-plotly .plot-container").wait_for(timeout=60000)
    page.wait_for_timeout(1200)
    capture(page, 4)
    with page.expect_download() as png_3d_download:
        plot_cards.nth(0).get_by_role("button", name="Export", exact=True).click()
        plot_cards.nth(0).get_by_role("button", name="PNG image", exact=True).click()
    assert png_3d_download.value.suggested_filename.endswith(".png")

    # Change only the map. The fitted FlowSOM object must be reused.
    map_select.select_option("diffusion-map")
    dialog.get_by_role("button", name="Run pipeline").click()
    dialog.get_by_text("FlowSOM → Diffusion map", exact=True).wait_for(timeout=120000)
    dialog.get_by_text("Cluster reused", exact=True).wait_for()
    capture(page, 5)

    dialog.locator(".analysis-result-tabs").get_by_role("button", name="Cell identities").click()
    identity_panel = dialog.locator(".analysis-identity-panel")
    identity_panel.wait_for()
    identity_cards = identity_panel.locator(".analysis-identity-card")
    assert identity_cards.count() == 2
    assert identity_cards.nth(0).get_by_role("button", name="CD3", exact=True).count() == 2
    assert identity_cards.nth(0).get_by_role("button", name="CD3-A", exact=True).count() == 0
    identity_panel.get_by_role("button", name="Cell identity settings").click()
    settings = identity_panel.get_by_role("dialog", name="Cell identity settings menu")
    settings.wait_for()
    assert settings.get_by_label("Minimum identity match").is_visible()
    assert settings.get_by_label("Minimum lead over second choice").is_visible()
    assert settings.get_by_label("Marker sensitivity").is_visible()
    assert "A cell must pass both confidence checks" in settings.inner_text()
    capture(page, 6)
    identity_panel.get_by_role("button", name="Cell identity settings").click()

    identity_panel.get_by_role("button", name="Load immune template (2)", exact=True).click()
    identity_names = [field.input_value() for field in identity_panel.get_by_label("Identity name").all()]
    assert identity_names == ["B cell", "NK cell"]
    identity_panel.get_by_role("button", name="Annotate cells", exact=True).click()
    identity_panel.get_by_text("assigned", exact=False).wait_for(timeout=60000)
    dialog.locator(".analysis-result-tabs").get_by_role("button", name="Plots", exact=True).click()
    assert "Predicted identity" in dialog.get_by_label("Result color").locator("option").all_inner_texts()
    capture(page, 7)

    # Trajectory methods declare prerequisites and cannot run before a root event is selected.
    dialog.get_by_role("button", name="Clear result", exact=True).click()
    assert dialog.locator(".analysis-result-tabs").get_by_role("button", name="Cell identities").is_disabled()
    dialog.get_by_role("button", name="Trajectory", exact=True).click()
    trajectory_select = dialog.get_by_label("Trajectory method")
    trajectory_select.wait_for()
    enabled_trajectories = {
        option.get_attribute("value")
        for option in trajectory_select.locator("option").all()
        if option.get_attribute("value") and option.is_enabled()
    }
    assert {"dpt", "slingshot", "tscan", "palantir", "paga-dpt", "wanderlust", "wishbone"}.issubset(enabled_trajectories)
    trajectory_select.select_option("slingshot")
    assert dialog.get_by_label("Trajectory clustering method").is_visible()
    assert dialog.get_by_label("Trajectory dimensional-reduction method").is_visible()
    assert dialog.get_by_role("button", name="Run trajectory").is_disabled()
    assert "root" in dialog.locator(".analysis-run-status").inner_text().lower()
    trajectory_select.select_option("wishbone")
    wishbone_settings = dialog.locator(".analysis-advanced-settings")
    wishbone_settings.locator("summary").click()
    assert wishbone_settings.get_by_text("Branch confidence", exact=True).is_visible()
    capture(page, 8)

    assert "[object Object]" not in page.locator("body").inner_text()
    browser.close()

if errors:
    raise AssertionError("Browser console errors:\n" + "\n".join(errors))
