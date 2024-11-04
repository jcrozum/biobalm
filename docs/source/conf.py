# type: ignore
# # Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = "biobalm"
copyright = "2024, Trinh, Park, Pastva, Rozum"
author = "Trinh, Park, Pastva, Rozum"
release = "0.4.2"

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "sphinx.ext.githubpages",
    "sphinx.ext.linkcode",
    "sphinx.ext.intersphinx",
    "sphinx.ext.napoleon",
    "sphinx_copybutton",
    "sphinx_toolbox.more_autodoc.genericalias",
    "autoapi.extension",
]
default_role = "code"
templates_path = ["_templates"]
exclude_patterns = []


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = "pydata_sphinx_theme"
html_static_path = ["_static"]
html_css_files = ["css/custom.css"]
html_theme_options = {
    "show_nav_level": 9,
    "navigation_with_keys": True,
    "show_toc_level": 9,
    "icon_links": [
        {
            # Label for this link
            "name": "GitHub",
            # URL where the link will redirect
            "url": "https://github.com/jcrozum/biobalm",  # required
            # Icon class (if "type": "fontawesome"), or path to local image (if "type": "local")
            "icon": "fa-brands fa-square-github",
            # The type of image to be used (see below for details)
            "type": "fontawesome",
        }
    ],
}

# -- Extension configuration -------------------------------------------------
autoapi_dirs = ["../../biobalm"]
autoapi_options = [
    "members",
    "undoc-members",
    "show-inheritance",
    "show-module-summary",
    # "special-members",
    "imported-members",
]
autoapi_python_class_content = "both"

intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "networkx": ("https://networkx.org/documentation/stable/", None),
}

# autodoc_type_aliases = {
#     "BooleanSpace": "BooleanSpace",
#     "SubspaceSuccession": "SubspaceSuccession",
#     "ControlOverrides": "ControlOverrides",
# }
autosummary_generate = True
autosummary_imported_members = True
autodoc_default_options = {
    "autosummary": True,
}
numpydoc_show_class_members = False
numpydoc_class_members_toctree = False
copybutton_prompt_text = r">>> |\.\.\. "
copybutton_prompt_is_regexp = True

automodapi_inheritance_diagram = False


def linkcode_resolve(domain, info):
    if domain != "py":
        return None
    if not info["module"]:
        return None

    if info["module"] == "biobalm":  # to handle top-level exports
        if info["fullname"].startswith("SuccessionDiagram"):
            info["module"] = "biobalm.succession_diagram"

    filename = info["module"].replace(".", "/")
    return f"https://github.com/jcrozum/biobalm/blob/main/{filename}.py"
