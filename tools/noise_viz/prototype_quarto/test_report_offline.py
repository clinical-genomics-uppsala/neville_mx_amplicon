"""PROTOTYPE - throwaway. Drive the Quarto prototype report offline in Chrome and check the settled interactions.

Usage (needs Playwright and a local Google Chrome):
    python test_report_offline.py <absolute path to noise_report.html> <folder for screenshots>

Every http(s) request is blocked, so the report must work with no network. Prints load and tab-draw times, the
plot state before and after each interaction, the blocked requests and console errors, and saves screenshots
q_*.png. It checks nothing itself: read the printout.
"""
import json
import sys
import time

from playwright.sync_api import sync_playwright

HTML, SHOTS = sys.argv[1], sys.argv[2]
blocked, console = [], []


def plot_state(page, name):
    """Read the current state of one amplicon's figure and variant table from the page.

    Args:
        page: Playwright page with the report open.
        name: amplicon name, e.g. "TP53_D2".

    Returns:
        None when the amplicon's figure is not drawn, else a dict:
            range: x-axis range; shapes: number of layout shapes (reference lines + highlight);
            highlight: whether a highlight rectangle is shown;
            rollNames: names of the rolling-mean traces (they include the window size);
            lowPoints: number of low-depth points shown;
            hp, hpName: number of homopolymers shown, and the trace name;
            varPoints: number of variant markers and labels shown;
            selRows: positions of the selected table rows; hiddenRows: number of hidden table rows.
    """
    return page.evaluate("""name => {
        const gd = document.getElementById(`nv-plot-${name}`);
        if (!gd || !gd.data) return null;
        const t = role => gd.data.filter(d => d.meta && d.meta.role === role);
        return {
            range: gd.layout.xaxis.range, shapes: gd.layout.shapes.length,
            highlight: gd.layout.shapes.some(s => s.type === "rect"),
            rollNames: t("roll").map(d => d.name),
            lowPoints: t("low").reduce((n, d) => n + d.y.filter(v => v != null).length, 0),
            hp: t("hp").map(d => d.x.length), hpName: t("hp").map(d => d.name),
            varPoints: t("var").reduce((n, d) => n + (d.x || []).length, 0),
            selRows: [...document.querySelectorAll(`#nv-amp-${name} tr.v.sel`)].map(r => r.dataset.pos),
            hiddenRows: [...document.querySelectorAll(`#nv-amp-${name} tr.v`)].filter(r => r.style.display === "none").length,
        };
    }""", name)


def drawn(page):
    """Return the element ids of the figures currently drawn (the report draws only the visible tab)."""
    return page.evaluate("[...document.querySelectorAll('.nv-plot')].filter(d => d.data).map(d => d.id)")


def set_input(page, label, value):
    """Set one report control as a user would, then wait 0.6 s for the figure to update.

    Args:
        page: Playwright page with the report open.
        label: start of the control's label text, e.g. "Rolling window".
        value: new value; True/False for a checkbox, a number or text otherwise.

    Finds the input inside the label, sets it and fires an "input" event, which nv.html listens to.
    """
    page.evaluate("""([label, value]) => {
        const lab = [...document.querySelectorAll('.nv-controls label')].find(l => l.textContent.trim().startsWith(label));
        const box = lab;
        const inp = box.querySelector('input[type=checkbox]') || box.querySelector('input[type=number]') || box.querySelector('input');
        if (inp.type === 'checkbox') inp.checked = value; else inp.value = value;
        inp.dispatchEvent(new Event('input', {bubbles: true}));
    }""", [label, value])
    page.wait_for_timeout(600)


with sync_playwright() as p:
    browser = p.chromium.launch(channel="chrome", headless=True)
    page = browser.new_page(viewport={"width": 1600, "height": 1300})
    page.route("http*://**/*", lambda r: (blocked.append(r.request.url), r.abort()))  # offline
    page.on("console", lambda m: console.append(f"{m.type}: {m.text}") if m.type in ("error", "warning") else None)
    page.on("pageerror", lambda e: console.append(f"pageerror: {e}"))
    t0 = time.time()
    page.goto(f"file://{HTML}")
    page.wait_for_function("document.querySelector('.nv-plot .main-svg') !== null", timeout=60000)
    print(f"load to first figure: {time.time() - t0:.1f} s; drawn: {drawn(page)}")
    page.wait_for_selector(".nv-controls form", timeout=30000)
    print("controls:", page.evaluate("[...document.querySelectorAll('.nv-controls label')].map(l => l.textContent.trim())"))
    page.screenshot(path=f"{SHOTS}/q_first.png")

    # tab switching: every amplicon once, check only the visible one is drawn and time each draw
    names = page.evaluate("[...document.querySelectorAll('.nv-amp')].map(d => d.dataset.amp)")
    times = []
    for n in names:
        t0 = time.time()
        page.get_by_role("tab", name=n, exact=True).click()
        page.wait_for_function(f"(() => {{ const g = document.getElementById('nv-plot-{n}'); return g && g.data; }})()")
        times.append(round(time.time() - t0, 2))
    print("tab draw times (s):", dict(zip(names, times)), "| drawn after loop:", drawn(page))

    # TP53_D2: click the first table row (known variant), expect zoom to +/-60 bp and a highlight
    page.get_by_role("tab", name="TP53_D2", exact=True).click()
    page.wait_for_timeout(500)
    row = page.locator("#nv-amp-TP53_D2 tr.v").first
    pos = int(row.get_attribute("data-pos"))
    row.click()
    page.wait_for_timeout(800)
    s = plot_state(page, "TP53_D2")
    print(f"row click {pos}: range {s['range']}, highlight {s['highlight']}, selected rows {s['selRows']}")
    page.screenshot(path=f"{SHOTS}/q_row_click.png")

    # plot click on a variant marker (emitted as Plotly would) highlights its row
    other = page.evaluate("""() => { const gd = document.getElementById('nv-plot-TP53_D2');
        const i = gd.data.findIndex(d => d.meta && d.meta.role === 'var' && d.x && d.x.length > 1);
        const t = gd.data[i]; gd.emit('plotly_click', {points: [{data: t, customdata: t.customdata[1], x: t.x[1]}]}); return t.x[1]; }""")
    page.wait_for_timeout(600)
    print(f"marker click {other}: selected rows {plot_state(page, 'TP53_D2')['selRows']}")

    # controls
    before = plot_state(page, "TP53_D2")
    set_input(page, "Rolling window", 50)
    set_input(page, "Minimum depth", 25000)
    set_input(page, "Min homopolymer", 8)
    set_input(page, "Reference lines", "2")
    after = plot_state(page, "TP53_D2")
    for k in ("rollNames", "lowPoints", "hp", "hpName", "shapes", "range"):
        print(f"  {k}: {before[k]} -> {after[k]}")
    page.screenshot(path=f"{SHOTS}/q_controls.png")

    page.get_by_role("tab", name="IDH1_1", exact=True).click()
    page.wait_for_timeout(600)
    b = plot_state(page, "IDH1_1")
    set_input(page, "Hide germline", True)
    a = plot_state(page, "IDH1_1")
    print(f"germline toggle IDH1_1: variant points {b['varPoints']} -> {a['varPoints']}, hidden rows {b['hiddenRows']} -> {a['hiddenRows']}",
          "| note:", page.inner_text("#nv-amp-IDH1_1 .nv-hidden"))
    set_input(page, "Germline VAF", 5)
    print("germline VAF >= 5%: variant points", plot_state(page, "IDH1_1")["varPoints"], "hidden rows", plot_state(page, "IDH1_1")["hiddenRows"])
    page.screenshot(path=f"{SHOTS}/q_germline.png")

    page.get_by_role("tab", name="Read quality", exact=True).click()
    page.wait_for_timeout(800)
    page.screenshot(path=f"{SHOTS}/q_readq.png")
    print("blocked network requests:", sorted(set(blocked)))
    print("console errors/warnings:", json.dumps(console[:10], indent=1))
    browser.close()
