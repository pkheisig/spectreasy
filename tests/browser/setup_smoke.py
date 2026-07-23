import sys
from pathlib import Path

from playwright.sync_api import sync_playwright


port = int(sys.argv[1])
capture = Path(sys.argv[2])
with sync_playwright() as playwright:
    browser = playwright.chromium.launch(headless=True)
    page = browser.new_page(viewport={"width": 1440, "height": 1050}, device_scale_factor=1)
    page.goto(f"http://127.0.0.1:{port}/?setup=1", wait_until="networkidle")
    page.get_by_role("heading", name="Setup", exact=True).wait_for()
    assert page.locator(".setup-panel").count() == 4
    assert page.get_by_role("heading", name="Install analysis methods", exact=True).is_visible()
    assert "spectreasy::install_analysis_dependencies()" in page.locator(".setup-analysis-panel").inner_text()
    assert "never changes system Python" in page.locator(".setup-analysis-panel").inner_text()
    assert "macOS app" not in page.locator("body").inner_text()
    assert "Windows app" not in page.locator("body").inner_text()
    capture.parent.mkdir(parents=True, exist_ok=True)
    page.screenshot(path=str(capture), full_page=True)
    browser.close()
