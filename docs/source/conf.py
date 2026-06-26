"""Sphinx configuration for the burnsidelib documentation."""

from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

# Some Sage extension classes expect the main Sage namespace to be initialized
# before they are imported by autodoc.
import sage.all  # noqa: E402,F401

project = "burnsidelib"
author = "burnsidelib contributors"
copyright = "2026, burnsidelib contributors"
release = "0.1"

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.mathjax",
    "sphinx_copybutton",
]

autodoc_member_order = "bysource"
autodoc_typehints = "description"

templates_path = ["_templates"]
exclude_patterns = []
highlight_language = "python"

html_theme = "furo"
html_title = "burnsidelib documentation"
html_static_path = ["_static"]
html_css_files = ["custom.css"]
html_theme_options = {
    "source_repository": "https://github.com/diracdeltafunk/burnsidelib/",
    "source_branch": "master",
    "source_directory": "docs/source/",
}

copybutton_prompt_text = r">>> |\.\.\. |sage: "
copybutton_prompt_is_regexp = True
