import tempfile
from pathlib import Path
from unittest import TestCase

from pling.html_visualisations import CYTOSCAPE_CDN_URL, repair_cytoscape_html_outputs


class TestHtmlVisualisations(TestCase):
    def test_adds_local_cytoscape_source_before_extensions(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            graph_dir = root / "visualisations" / "communities" / "graphs"
            libs_dir = root / "visualisations" / "communities" / "libs"
            graph_dir.mkdir(parents=True)
            libs_dir.mkdir(parents=True)
            (libs_dir / "cytoscape.min.js").write_text("window.cytoscape = function () {};")
            html_path = graph_dir / "community_0.html"
            html_path.write_text(
                "\n".join(
                    [
                        "<html>",
                        "<head>",
                        '    <script src="../libs/jquery-3.6.0.min.js"></script>',
                        '    <script src="../libs/cytoscape-euler.js"></script>',
                        "</head>",
                        "<body>",
                        "    <script>cytoscape({});</script>",
                        "</body>",
                        "</html>",
                    ]
                )
            )

            self.assertEqual(repair_cytoscape_html_outputs(root), 1)
            repaired = html_path.read_text()
            self.assertIn('<script src="../libs/cytoscape.min.js"></script>', repaired)
            self.assertLess(
                repaired.index("cytoscape.min.js"),
                repaired.index("cytoscape-euler.js"),
            )
            self.assertEqual(repair_cytoscape_html_outputs(root), 0)

    def test_uses_cdn_when_no_local_libs_are_available(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            html_path = root / "graph.html"
            html_path.write_text(
                "<html><head></head><body><script>cytoscape({});</script></body></html>"
            )

            self.assertEqual(repair_cytoscape_html_outputs(root), 1)
            self.assertIn(CYTOSCAPE_CDN_URL, html_path.read_text())

    def test_leaves_non_cytoscape_html_unchanged(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            html_path = root / "index.html"
            html_path.write_text("<html><body>No graph here.</body></html>")

            self.assertEqual(repair_cytoscape_html_outputs(root), 0)
            self.assertEqual(html_path.read_text(), "<html><body>No graph here.</body></html>")
