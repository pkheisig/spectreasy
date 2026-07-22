import sys
from pathlib import Path
from playwright.sync_api import sync_playwright

port = int(sys.argv[1])
output = Path(sys.argv[2])
errors = []

with sync_playwright() as playwright:
    browser = playwright.chromium.launch(headless=True)
    page = browser.new_page(viewport={"width": 1600, "height": 1000}, device_scale_factor=1)
    page.on("console", lambda message: errors.append(message.text) if message.type == "error" else None)
    page.on("pageerror", lambda error: errors.append(str(error)))
    page.goto(
        f"http://127.0.0.1:{port}/?api=http%3A%2F%2F127.0.0.1%3A{port}#token=analysis-v2-smoke-token",
        wait_until="networkidle",
    )
    page.locator("body").wait_for()
    page.locator(".rail-subitem", has_text="Samples").click()
    page.get_by_role("button", name="Analyze samples").click()
    page.get_by_text("v2 workspace", exact=True).wait_for()
    page.get_by_role("button", name="Rectangle").click()
    canvas = page.locator(".analysis-plot-surface canvas").first
    box = canvas.bounding_box()
    assert box is not None
    page.mouse.move(box["x"] + box["width"] * 0.30, box["y"] + box["height"] * 0.30)
    page.mouse.down()
    page.mouse.move(box["x"] + box["width"] * 0.70, box["y"] + box["height"] * 0.72)
    page.mouse.up()
    dialog = page.locator(".analysis-gate-dialog")
    dialog.get_by_label("Name").fill("Browser gate")
    dialog.get_by_label("Role").select_option("negative")
    dialog.get_by_role("button", name="Create population").click()
    page.locator(".analysis-population-row", has_text="Browser gate").first.wait_for()
    page.get_by_role("button", name="Add plot").click()
    assert page.locator(".analysis-plot-card").count() == 2
    assert page.locator(".analysis-populations").is_visible()
    assert page.get_by_role("button", name="Discover", exact=True).is_visible()
    page.locator(".analysis-save-state.is-saved").wait_for()
    assert "[object Object]" not in page.locator("body").inner_text()
    page.screenshot(path=str(output), full_page=True)
    browser.close()

if errors:
    raise AssertionError("Browser console errors:\n" + "\n".join(errors))
