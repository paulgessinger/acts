#!/usr/bin/env python3

from __future__ import annotations

import argparse
import ast
import keyword
from pathlib import Path
import re
import sys


TYPEVAR_LINE = '_ReadDataHandleT = typing.TypeVar("_ReadDataHandleT")'

READ_DATA_HANDLE_BLOCK = [
    "class ReadDataHandle(typing.Generic[_ReadDataHandleT]):",
    "    def __call__(self, whiteboard: WhiteBoard) -> _ReadDataHandleT:",
    "        ...",
    "    def __init__(self, parent: typing.Any, type: type[_ReadDataHandleT], name: str) -> None:",
    "        ...",
    "    def initialize(self, key: str) -> None:",
    "        ...",
]


_INIT_OVERLOAD_RE = re.compile(
    r"__init__\(\s*self:\s*[^,]+,\s*config:\s*[^,\)]+(?P<rest>.*)\)\s*->\s*None"
)
_TYPED_PARAM_RE = re.compile(
    r"(?P<name>[A-Za-z_][A-Za-z0-9_]*)\s*:\s*(?P<ann>.+?)(?:\s*=\s*(?P<default>.+))?$"
)
_ANN_ASSIGN_LINE_RE = re.compile(
    r"^(?P<indent>\s*)(?P<name>[A-Za-z_][A-Za-z0-9_]*)(?P<sep>\s*:\s*)(?P<rest>.*)$"
)
_DEF_LINE_RE = re.compile(
    r"^(?P<indent>\s*)def\s+(?P<name>[A-Za-z_][A-Za-z0-9_]*)\s*\(.*$"
)


def _sanitize_for_ast_parse(text: str) -> str:
    """Make generated stubs parseable for AST inspection.

    pybind11-stubgen can emit class var annotations with keyword targets like
    `None: typing.ClassVar[...]`, which are invalid syntax. We rewrite these to
    `None_: ...` only in an in-memory copy used for AST parsing.
    """

    lines = []
    for line in text.splitlines():
        dm = _DEF_LINE_RE.match(line)
        if dm is not None and "::" in line:
            lines.append(f"{dm.group('indent')}def {dm.group('name')}(*args, **kwargs):")
            continue

        m = _ANN_ASSIGN_LINE_RE.match(line)
        if m is None:
            lines.append(line)
            continue

        name = m.group("name")
        if keyword.iskeyword(name):
            lines.append(
                f"{m.group('indent')}{name}_{m.group('sep')}{m.group('rest')}"
            )
        else:
            lines.append(line)

    return "\n".join(lines) + "\n"


def _parse_stub_ast(text: str) -> ast.Module:
    try:
        return ast.parse(text)
    except SyntaxError:
        return ast.parse(_sanitize_for_ast_parse(text))


def _split_params(rest: str) -> list[str]:
    params = []
    depth = 0
    current = []
    for ch in rest:
        if ch in "([{":
            depth += 1
        elif ch in ")]}":
            depth = max(0, depth - 1)
        elif ch == "," and depth == 0:
            token = "".join(current).strip()
            if token:
                params.append(token)
            current = []
            continue
        current.append(ch)
    token = "".join(current).strip()
    if token:
        params.append(token)
    return params


def _is_valid_annotation(ann: str) -> bool:
    try:
        ast.parse(f"x: {ann}\n")
    except SyntaxError:
        return False
    return True


def _is_valid_expression(expr: str) -> bool:
    try:
        ast.parse(expr, mode="eval")
    except SyntaxError:
        return False
    return True


def _sanitize_annotation(ann: str) -> str:
    return ann if _is_valid_annotation(ann) else "typing.Any"


def _sanitize_default(default: str | None) -> str | None:
    if default is None:
        return None
    return default if _is_valid_expression(default) else "..."


def _parse_ctor_overload_doc(doc: str) -> list[list[tuple[str, str, str | None]]]:
    overloads: list[list[tuple[str, str, str | None]]] = []
    for raw_line in doc.splitlines():
        line = raw_line.strip()
        if "__init__(" not in line or "config:" not in line:
            continue

        m = _INIT_OVERLOAD_RE.search(line)
        if m is None:
            continue

        rest = m.group("rest").strip()
        if rest.startswith(","):
            rest = rest[1:].strip()

        params: list[tuple[str, str, str | None]] = []
        for token in _split_params(rest):
            pm = _TYPED_PARAM_RE.fullmatch(token)
            if pm is None:
                continue
            ann = _sanitize_annotation(pm.group("ann").strip())
            default = _sanitize_default(pm.group("default"))
            params.append(
                (pm.group("name"), ann, default)
            )

        if params not in overloads:
            overloads.append(params)

    return overloads


def _is_property_decorator(node: ast.expr) -> bool:
    return isinstance(node, ast.Name) and node.id == "property"


def _collect_config_fields(cfg_cls: ast.ClassDef, text: str) -> list[tuple[str, str]]:
    fields: list[tuple[str, str]] = []
    seen: set[str] = set()

    for stmt in cfg_cls.body:
        if (
            isinstance(stmt, ast.AnnAssign)
            and isinstance(stmt.target, ast.Name)
            and stmt.annotation is not None
        ):
            name = stmt.target.id
            ann = ast.get_source_segment(text, stmt.annotation)
            if ann is None or not name.isidentifier() or name in seen:
                continue
            fields.append((name, ann))
            seen.add(name)

    for stmt in cfg_cls.body:
        if (
            isinstance(stmt, ast.FunctionDef)
            and any(_is_property_decorator(d) for d in stmt.decorator_list)
            and stmt.returns is not None
        ):
            name = stmt.name
            ann = ast.get_source_segment(text, stmt.returns)
            if ann is None or not name.isidentifier() or name in seen:
                continue
            fields.append((name, ann))
            seen.add(name)

    return fields


def _render_params(params: list[tuple[str, str, str | None]]) -> str:
    rendered = []
    for name, ann, default in params:
        if default is None:
            rendered.append(f"{name}: {ann}")
        else:
            rendered.append(f"{name}: {ann} = {default}")
    return ", ".join(rendered)


def _render_rewritten_init(
    config_fields: list[tuple[str, str]],
    ctor_overloads: list[list[tuple[str, str, str | None]]],
) -> list[str]:
    block: list[str] = []

    if ctor_overloads:
        for params in ctor_overloads:
            config_tail = _render_params(params)
            if config_tail:
                block.append("    @typing.overload")
                block.append(
                    f"    def __init__(self, config: Config, {config_tail}) -> None:"
                )
                block.append("        ...")
            else:
                block.append("    @typing.overload")
                block.append("    def __init__(self, config: Config) -> None:")
                block.append("        ...")

    if config_fields:
        kwargs_fields = [f"{name}: {ann} = ..." for name, ann in config_fields]
        for params in ctor_overloads or [[]]:
            required_ctor = [
                f"{name}: {ann}" for name, ann, default in params if default is None
            ]
            optional_ctor = [
                f"{name}: {ann} = {default}"
                for name, ann, default in params
                if default is not None
            ]
            kwargs = required_ctor + kwargs_fields + optional_ctor
            block.append("    @typing.overload")
            block.append(f"    def __init__(self, *, {', '.join(kwargs)}) -> None:")
            block.append("        ...")

    block.append("    def __init__(self, *args, **kwargs):")
    block.append("        ...")
    return block


def _patch_config_adapted_constructors(text: str) -> str:
    tree = _parse_stub_ast(text)
    lines = text.splitlines()
    replacements: list[tuple[int, int, list[str]]] = []

    for node in tree.body:
        if not isinstance(node, ast.ClassDef):
            continue

        config_cls = None
        init_fns: list[ast.FunctionDef] = []
        for stmt in node.body:
            if isinstance(stmt, ast.ClassDef) and stmt.name == "Config":
                config_cls = stmt
            elif isinstance(stmt, ast.FunctionDef) and stmt.name == "__init__":
                init_fns.append(stmt)

        if config_cls is None or not init_fns:
            continue

        init_fn = init_fns[-1]

        if (
            len(init_fn.args.args) != 1
            or init_fn.args.vararg is None
            or init_fn.args.kwarg is None
        ):
            continue

        cfg_fields = _collect_config_fields(config_cls, text)
        ctor_doc = ast.get_docstring(init_fn, clean=False) or ""
        ctor_overloads = _parse_ctor_overload_doc(ctor_doc)

        if not cfg_fields and not ctor_overloads:
            continue

        block = _render_rewritten_init(cfg_fields, ctor_overloads)
        replacements.append(
            (init_fns[0].lineno - 1, init_fns[-1].end_lineno - 1, block)
        )

    for start, end, block in sorted(
        replacements, key=lambda item: item[0], reverse=True
    ):
        lines[start : end + 1] = block

    return "\n".join(lines) + "\n"


def _find_import_typing_lineno(tree: ast.Module) -> int | None:
    for node in tree.body:
        if isinstance(node, ast.Import):
            for alias in node.names:
                if alias.name == "typing" and alias.asname is None:
                    return node.lineno
    return None


def _find_class(tree: ast.Module, name: str) -> ast.ClassDef | None:
    for node in tree.body:
        if isinstance(node, ast.ClassDef) and node.name == name:
            return node
    return None


def _patch_read_data_handle(text: str) -> str:
    lines = text.splitlines()
    tree = _parse_stub_ast(text)

    read_handle_cls = _find_class(tree, "ReadDataHandle")
    if read_handle_cls is None:
        return text

    if read_handle_cls.end_lineno is None:
        return text

    if not any(
        isinstance(stmt, ast.FunctionDef) and stmt.name == "__call__"
        for stmt in read_handle_cls.body
    ):
        return text

    typing_import_lineno = _find_import_typing_lineno(tree)
    if typing_import_lineno is None:
        return text

    if TYPEVAR_LINE not in lines:
        insert_at = typing_import_lineno
        lines.insert(insert_at, TYPEVAR_LINE)

        # Class line numbers are 1-based and shift by one when insertion
        # happens above.
        if read_handle_cls.lineno > insert_at:
            class_start = read_handle_cls.lineno
            class_end = read_handle_cls.end_lineno
        else:
            class_start = read_handle_cls.lineno - 1
            class_end = read_handle_cls.end_lineno - 1
    else:
        class_start = read_handle_cls.lineno - 1
        class_end = read_handle_cls.end_lineno - 1

    lines[class_start : class_end + 1] = READ_DATA_HANDLE_BLOCK
    return "\n".join(lines) + "\n"


def patch_stub_file(path: Path) -> bool:
    text = path.read_text()
    updated = _patch_read_data_handle(text)
    updated = _patch_config_adapted_constructors(updated)

    if updated == text:
        return False

    path.write_text(updated)
    return True


def _iter_stub_paths(root: Path, explicit_stubs: list[Path]) -> list[Path]:
    seen: set[Path] = set()
    paths: list[Path] = []

    for stub in explicit_stubs:
        resolved = stub.resolve()
        if resolved in seen or not stub.exists():
            continue
        seen.add(resolved)
        paths.append(stub)

    if root.exists():
        for stub in sorted(root.rglob("*.pyi")):
            resolved = stub.resolve()
            if resolved in seen:
                continue
            seen.add(resolved)
            paths.append(stub)

    return paths


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Patch generated pybind11 stubs for improved typing."
    )
    parser.add_argument(
        "--root",
        type=Path,
        default=Path("acts"),
        help="Root directory under which .pyi files are patched recursively.",
    )
    parser.add_argument(
        "--stub",
        type=Path,
        action="append",
        default=[],
        help="Additional explicit .pyi files to patch.",
    )
    args = parser.parse_args()

    for stub in _iter_stub_paths(args.root, args.stub):
        try:
            patch_stub_file(stub)
        except SyntaxError as exc:
            print(f"{stub}: syntax error while patching stubs: {exc}", file=sys.stderr)


if __name__ == "__main__":
    main()
