import sys
from pathlib import Path

from playwright.sync_api import sync_playwright


port = int(sys.argv[1])
captures = [Path(value) for value in sys.argv[2:]]
errors: list[str] = []


def capture(page, index: int) -> None:
    if index < len(captures):
        captures[index].parent.mkdir(parents=True, exist_ok=True)
        page.screenshot(path=str(captures[index]), full_page=False)


with sync_playwright() as playwright:
    browser = playwright.chromium.launch(headless=True)
    page = browser.new_page(viewport={"width": 820, "height": 900}, device_scale_factor=1)
    page.on("console", lambda message: errors.append(message.text) if message.type == "error" else None)
    page.on("pageerror", lambda error: errors.append(str(error)))
    page.goto(
        f"http://127.0.0.1:{port}/?api=http%3A%2F%2F127.0.0.1%3A{port}#token=analysis-v2-smoke-token",
        wait_until="networkidle",
    )
    prompt = page.locator(".project-initialize-confirm")
    if prompt.count() and prompt.is_visible():
        prompt.locator(".button-ghost", has_text="Not now").click()
    page.get_by_role("button", name="Samples").click()
    assert page.get_by_role("button", name="Analyze samples").count() == 0
    page.get_by_role("button", name="Other tools").click()
    page.locator(".rail-subitem", has_text="Population analysis").click()
    page.get_by_role("button", name="Analyze population", exact=True).wait_for(timeout=60000)
    workspace = page.locator(".analysis-shell")
    workspace.wait_for()
    workspace_box = workspace.bounding_box()
    assert workspace_box["x"] == 0 and workspace_box["width"] <= 820
    capture(page, 0)
    assert page.evaluate("window.scrollTo(500, 0); window.scrollX === 0")

    page.get_by_role("button", name="Analyze population", exact=True).click()
    dialog = page.locator(".analysis-method-dialog")
    dialog.wait_for()
    assert dialog.bounding_box()["width"] <= 820
    cluster = dialog.locator(".analysis-pipeline-stage").nth(0)
    cluster.get_by_role("checkbox", name="Skip").check()
    dialog.get_by_label("Dimensional-reduction method").select_option("pca")
    dialog.get_by_role("button", name="Run pipeline").click()
    dialog.get_by_text("PCA", exact=True).last.wait_for(timeout=120000)
    plot = dialog.locator(".analysis-result-plot")
    plot.wait_for()
    box = plot.bounding_box()
    assert box is not None and abs(box["width"] - box["height"]) < 2
    assert box["width"] <= 380
    assert page.evaluate("window.scrollTo(500, 0); window.scrollX === 0")
    capture(page, 1)
    browser.close()

if errors:
    raise AssertionError("Browser console errors:\n" + "\n".join(errors))
