def process_component_docstring(app, what, name, obj, options, lines):
    from cratermaker.components.counting import Counting
    from cratermaker.components.crater import Crater
    from cratermaker.components.morphology import Morphology
    from cratermaker.components.production import Production
    from cratermaker.components.projectile import Projectile
    from cratermaker.components.scaling import Scaling
    from cratermaker.components.surface import Surface
    from cratermaker.components.target import Target

    component_classes = {
        "Scaling": Scaling,
        "Production": Production,
        "Morphology": Morphology,
        "Surface": Surface,
        "Projectile": Projectile,
        "Target": Target,
        "Counting": Counting,
        "Crater": Crater,
    }

    if name.endswith(".maker"):
        parent_name = name.rsplit(".", 1)[0]
        cls_name = parent_name.split(".")[-1].capitalize()
        cls = component_classes.get(cls_name)
        if cls:
            try:
                available = cls.available()
                if available:
                    lines.append("")
                    lines.append("**Available models:**")
                    lines.append("")
                    for item in available:
                        lines.append(f"- ``{item}``")
            except Exception as e:
                lines.append("")
                lines.append(f"*Error listing available components: {e}*")


def setup(app):
    app.connect("autodoc-process-docstring", process_component_docstring)
