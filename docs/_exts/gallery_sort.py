"""Sort key for the examples of the sphinx-gallery (docs/conf.py, sphinx_gallery_conf["within_subsection_order"])."""

import re


def example_number(filename: str) -> tuple:
    """
    Order gallery examples by their numeric prefix, so that 1.2 comes before 1.10.

    Parameters
    ----------
    filename : str
        The file name of the example, beginning with "<section>.<number>" (e.g. "1.10-save_actions.py").

    Returns
    -------
    tuple
        The section and example numbers, then the file name; files without a numeric prefix sort after the numbered ones.
    """
    match = re.match(r"(\d+)\.(\d+)", filename)
    if match is None:
        return (float("inf"), float("inf"), filename)
    return (int(match.group(1)), int(match.group(2)), filename)
