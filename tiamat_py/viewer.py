from __future__ import annotations

from dataclasses import dataclass, field
import json
from math import cos, pi, radians, sin, sqrt, tan
from pathlib import Path
import subprocess
import sys
import tempfile
import tkinter as tk
from tkinter import colorchooser, filedialog, messagebox, ttk

from .io_dna import load_dna_project
from .io_json import load_json_project, save_json_project
from .io_pdb import load_pdb_sugar_center_project
from .model import BASE_COLORS, Base, TiamatProject, basepair_state


VIEW_LAYOUT = (
    ("xz", "XZ", "xz", 0, 0),
    ("view3d", "3D", "3d", 0, 1),
    ("yz", "YZ", "yz", 1, 0),
    ("xy", "XY", "xy", 1, 1),
)

VIEW_BADGES = {
    "xz": "T",
    "yz": "S",
    "xy": "F",
}

RENDER_VIEW_OPTIONS = (
    ("Perspective", "view3d"),
    ("Top", "xz"),
    ("Front", "xy"),
    ("Side", "yz"),
)

RENDER_DESTINATION_OPTIONS = ("File", "Clipboard")

CONSTRAINT_COLUMNS = (
    ("b_median", "B-Type Median"),
    ("a_median", "A-Type Median"),
)

CUSTOM_CONSTRAINT_COLUMNS = (
    ("x_median", "X-Type Median"),
)

CONSTRAINT_ROWS = (
    ("rotation_per_base_pair_deg", "Rotation per Base Pair (deg)"),
    ("helicity_bp_per_turn", "Helicity"),
    ("rise_per_base_pair_nm", "Rise per Base Pair (nm)"),
    ("diameter_nm", "Diameter (nm)"),
    ("inclination_deg", "Inclination (deg)"),
    ("minor_groove_angle_deg", "Minor Groove Angle (deg)"),
)

DEFAULT_CONSTRAINT_VALUES = {
    "rotation_per_base_pair_deg": {
        "b_median": "-34.2857",
        "a_median": "-32.7",
    },
    "helicity_bp_per_turn": {
        "b_median": "10.50",
        "a_median": "11.01",
    },
    "rise_per_base_pair_nm": {
        "b_median": "0.332",
        "a_median": "0.29",
    },
    "diameter_nm": {
        "b_median": "2",
        "a_median": "2.3",
    },
    "inclination_deg": {
        "b_median": "0",
        "a_median": "19",
    },
    "minor_groove_angle_deg": {
        "b_median": "135.92",
        "a_median": "106.17",
    },
}

AXIS_COLORS = {
    "x": "#ef4444",
    "y": "#22c55e",
    "z": "#60a5fa",
}

COLOR_SCHEMES = {
    "dark": {
        "canvas_bg": "#101417",
        "viewport_bg": "#0f172a",
        "canvas_border": "#334155",
        "sash_bg": "#475569",
        "grid_line": "#334155",
        "selection_outline": "#ffffff",
        "base_outline": "#111827",
        "label_fg": "#f8fafc",
        "badge_fg": "#e5e7eb",
        "selection_box": "#93c5fd",
    },
    "bright": {
        "canvas_bg": "#ffffff",
        "viewport_bg": "#f8fafc",
        "canvas_border": "#cbd5e1",
        "sash_bg": "#cbd5e1",
        "grid_line": "#d1d5db",
        "selection_outline": "#111111",
        "base_outline": "#111111",
        "label_fg": "#111111",
        "badge_fg": "#111111",
        "selection_box": "#111111",
    },
}

SELECTION_MODES = (
    ("helix", "Select Helix", "select_helix", "1"),
    ("strand", "Select Strand", "select_strand", "2"),
    ("base_pair", "Select Base Pair", "select_base_pair", "3"),
    ("base", "Select Base", "select_base", "4"),
    ("box", "Select Box", "select_box", "5"),
)
CREATE_MODES = (
    ("create_helix", "Create Helix", "create_helix", "6"),
    ("create_base", "Create Free Strand", "create_base", "7"),
)
TRANSFORM_MODES = (
    ("move", "Move", "move", "8"),
    ("rotate", "Rotate", "rotate", "9"),
)
INTERACTION_MODES = SELECTION_MODES + CREATE_MODES + TRANSFORM_MODES
INTERACTION_MODE_LABELS = {mode: label for mode, label, _icon, _key in INTERACTION_MODES}
INTERACTION_MODE_ICONS = {mode: icon for mode, _label, icon, _key in INTERACTION_MODES}
TOOLBAR_MODE_ICONS = {icon for _mode, _label, icon, _key in INTERACTION_MODES}
MODE_KEYS = {shortcut: mode for mode, _label, _icon, shortcut in INTERACTION_MODES if shortcut}
TRANSFORM_MODE_NAMES = {mode for mode, _label, _icon, _key in TRANSFORM_MODES}
EMPTY_CLICK_CLEAR_MODES = {"helix", "strand", "base_pair", "base"}

TOOLBAR_HINTS = {
    "new_file": "New",
    "open_file": "Open",
    "save": "Save",
    "camera": "Render",
    "cut": "Cut",
    "copy": "Copy",
    "paste": "Paste",
    "delete": "Delete",
    "create_complementary": "Create Complementary",
    "delete_complementary": "Delete Complementary",
    "create_down": "Create Down",
    "delete_down": "Delete Down",
    "create_sticky_end": "Create Sticky End",
    "delete_sticky_end": "Delete Sticky End",
    "change_backbone_color": "Change Backbone Color",
    "reset_backbone_color": "Reset Backbone Color",
    "set_sequence": "Set Sequence",
    "reset_sequence": "Reset Sequence",
    "generate_sequence": "Generate Sequence",
    "undo": "Undo",
    "redo": "Redo",
    "strand_number": "Strand Number",
    "show_sticky_end": "Show Sticky Ends",
    "grid": "Grid",
    "select_helix": "Select Helix",
    "select_strand": "Select Strand",
    "select_base_pair": "Select Base Pair",
    "select_base": "Select Base",
    "select_box": "Select Box",
    "create_helix": "Create Helix",
    "create_base": "Create Free Strand",
    "move": "Move",
    "move_fixed": "Move Fixed",
    "rotate": "Rotate",
    "rotate_fixed": "Rotate Fixed",
}

LEGEND_ITEMS = (
    ("N", "#9ca3af"),
    ("A", "#22c55e"),
    ("T", "#ef4444"),
    ("U", "#a855f7"),
    ("C", "#3b82f6"),
    ("G", "#facc15"),
)

CREATE_MOLECULE_OPTIONS = ("DNA-DNA", "RNA-RNA", "DNA-RNA", "RNA-DNA")
CREATE_HELIX_FAMILY_OPTIONS = ("B-type", "A-type")

CREATE_STRUCTURE_TYPES = ("Helix", "Strand")

DEFAULT_CREATE_BASE_COUNT = 60

ABSOLUTE_ZOOM_MIN = 0.02
ABSOLUTE_ZOOM_MAX = 600.0
MAX_HISTORY_STEPS = 50
CLIPBOARD_RENDER_MAX_DIMENSION = 2400.0
IS_WINDOWS = sys.platform.startswith("win")
IS_MAC = sys.platform == "darwin"
ROTATE_DRAG_BUTTON = 3 if sys.platform == "darwin" else 2
PAN_DRAG_BUTTON = 2 if sys.platform == "darwin" else 3
VIEWER_STATE_PATH = Path.home() / ".tiamat_py_viewer_state.json"
OPENABLE_SUFFIXES = {".dna", ".json", ".dnajson", ".pdb"}


@dataclass
class Camera:
    rot_x: float = 20.0
    rot_y: float = -25.0
    rot_z: float = 0.0
    zoom: float = 0.0
    pan_x: float = 0.0
    pan_y: float = 0.0
    orientation: tuple[tuple[float, float, float], tuple[float, float, float], tuple[float, float, float]] | None = None


@dataclass
class ViewportState:
    key: str
    title: str
    mode: str
    canvas: tk.Canvas
    camera: Camera = field(default_factory=Camera)
    projected: dict[int, tuple[float, float, float]] = field(default_factory=dict)
    drag_anchor: tuple[int, int] | None = None
    drag_camera: Camera | None = None
    drag_button: int | None = None
    drag_moved: bool = False
    drag_positions: dict[int, tuple[float, float, float]] | None = None
    drag_center: tuple[float, float, float] | None = None
    drag_projected_center: tuple[float, float] | None = None
    drag_snapshot: dict[str, object] | None = None
    selection_box_start: tuple[int, int] | None = None
    selection_box_end: tuple[int, int] | None = None
    needs_fit: bool = True
    fit_zoom: float = 0.0
    scene_center: tuple[float, float, float] | None = None


@dataclass(frozen=True)
class FreeStrandEndpoint:
    point: tuple[float, float, float]
    base_index: int | None = None
    terminal_roles: frozenset[str] = frozenset()
    molecule: str | None = None


def _resource_root() -> Path:
    bundle_root = getattr(sys, "_MEIPASS", None)
    if bundle_root:
        return Path(bundle_root)
    return Path(__file__).resolve().parent.parent


def _shortcut_label(name: str) -> str:
    mac = {
        "new": "Cmd+N",
        "open": "Cmd+O",
        "save": "Cmd+S",
        "save_as": "Shift+Cmd+S",
        "cut": "Cmd+X",
        "copy": "Cmd+C",
        "paste": "Cmd+V",
        "undo": "Cmd+Z",
        "redo": "Shift+Cmd+Z",
        "select_all": "Cmd+A",
    }
    other = {
        "new": "Ctrl+N",
        "open": "Ctrl+O",
        "save": "Ctrl+S",
        "save_as": "Ctrl+Shift+S",
        "cut": "Ctrl+X",
        "copy": "Ctrl+C",
        "paste": "Ctrl+V",
        "undo": "Ctrl+Z",
        "redo": "Ctrl+Shift+Z",
        "select_all": "Ctrl+A",
    }
    labels = mac if IS_MAC else other
    return labels[name]


def _mark_viewport_for_fit(state: ViewportState, recenter: bool = False) -> None:
    state.needs_fit = True
    state.fit_zoom = 0.0
    state.scene_center = None
    if recenter:
        state.camera.pan_x = 0.0
        state.camera.pan_y = 0.0


class TiamatViewer:
    def __init__(self, root: tk.Tk, project: TiamatProject | None = None) -> None:
        self.root = root
        self.project = project
        self.selected_indices: set[int] = set()
        self.primary_selected: int | None = None
        self.viewports: dict[str, ViewportState] = {}
        self.menu_labels: dict[str, tk.Label] = {}
        self.menu_factories: dict[str, object] = {}
        self.native_menu_bar: tk.Menu | None = None
        self.native_menu_vars: list[tk.BooleanVar] = []
        self.menu_popup: tk.Frame | None = None
        self.submenu_popup: tk.Frame | None = None
        self._active_menu_name: str | None = None
        self._active_submenu_widgets: tuple[tk.Widget, ...] | None = None
        self.toolbar_icons: dict[str, tk.PhotoImage] = {}
        self.toolbar_disabled_icons: dict[str, tk.PhotoImage] = {}
        self.toolbar_buttons: dict[str, tk.Button] = {}
        self.toolbar_slots: dict[str, tk.Frame] = {}
        self.toolbar_hint: tk.Label | None = None
        self._sash_thickness = 8
        self._col_ratio = 0.5
        self._row_ratio = 0.5
        self._active_sash: str | None = None
        self.display_sphere_size = 4.0
        self.display_backbone_width = 2
        self.display_basepair_width = 1
        self.grid_extents = {
            "x_min": -5.0,
            "x_max": 5.0,
            "y_min": -5.0,
            "y_max": 5.0,
            "z_min": -5.0,
            "z_max": 5.0,
            "spacing": 1.0,
        }
        self.color_scheme = "dark"
        self.selection_mode = "base"
        self.pending_free_strand_start: FreeStrandEndpoint | None = None
        self.pending_free_strand_hover_view: str | None = None
        self.pending_free_strand_hover_canvas: tuple[int, int] | None = None
        self.show_strand_labels = True
        self.show_sticky_end = True
        self.show_grid = True
        self.constraints_values = _default_constraints_values()
        self.custom_constraints_values = _blank_custom_constraints_values()
        self.recent_files = self._load_recent_files()
        self.undo_stack: list[dict[str, object]] = []
        self.redo_stack: list[dict[str, object]] = []
        self.structure_clipboard: dict[str, object] | None = None
        self.structure_clipboard_paste_count = 0

        root.title("TiamatPy")
        root.geometry("1380x920")
        self._apply_window_icon()

        self._build_ui()
        self._apply_color_scheme()
        self._bind_events()
        self._reset_viewports()
        self._refresh_all()

    def _apply_window_icon(self) -> None:
        png_path = _resource_root() / "Toolbar_icons" / "exe_logo.png"
        ico_path = _resource_root() / "Toolbar_icons" / "exe_logo.ico"
        try:
            if png_path.exists():
                icon_image = tk.PhotoImage(file=str(png_path))
                self.root.iconphoto(True, icon_image)
                self._window_icon_image = icon_image
            if IS_WINDOWS and ico_path.exists():
                self.root.iconbitmap(default=str(ico_path))
        except tk.TclError:
            pass

    def _build_ui(self) -> None:
        self.main = ttk.Frame(self.root, padding=8)
        self.main.pack(fill=tk.BOTH, expand=True)

        self._initialize_menu_factories()
        if IS_WINDOWS:
            self._build_native_menu_bar()
        else:
            self._build_menu_strip()
        self._build_toolbar()

        self.viewport_host = tk.Frame(self.main, background="#0f172a", highlightthickness=0)
        self.viewport_host.pack(fill=tk.BOTH, expand=True, pady=(6, 0))
        self.viewport_host.bind("<Configure>", lambda _event: self._layout_viewports())

        self.vertical_sash = tk.Frame(
            self.viewport_host,
            background="#475569",
            cursor="sb_h_double_arrow",
            highlightthickness=0,
        )
        self.horizontal_sash = tk.Frame(
            self.viewport_host,
            background="#475569",
            cursor="sb_v_double_arrow",
            highlightthickness=0,
        )
        self.vertical_sash.bind("<ButtonPress-1>", lambda _event: self._active_sash_set("vertical"))
        self.horizontal_sash.bind("<ButtonPress-1>", lambda _event: self._active_sash_set("horizontal"))
        self.vertical_sash.bind("<B1-Motion>", self._drag_vertical_sash)
        self.horizontal_sash.bind("<B1-Motion>", self._drag_horizontal_sash)
        self.vertical_sash.bind("<ButtonRelease-1>", lambda _event: self._active_sash_set(None))
        self.horizontal_sash.bind("<ButtonRelease-1>", lambda _event: self._active_sash_set(None))

        for key, title, mode, row, column in VIEW_LAYOUT:
            container = tk.Frame(self.viewport_host, background="#0f172a", highlightthickness=0)

            canvas = tk.Canvas(
                container,
                background="#101417",
                highlightthickness=1,
                highlightbackground="#334155",
                highlightcolor="#93c5fd",
            )
            canvas.pack(fill=tk.BOTH, expand=True)

            self.viewports[key] = ViewportState(
                key=key,
                title=title,
                mode=mode,
                canvas=canvas,
                camera=_default_camera(mode),
            )

        self.status_var = tk.StringVar(value="No project loaded.")
        ttk.Label(self.main, textvariable=self.status_var, anchor="w").pack(fill=tk.X, pady=(6, 0))

    def _initialize_menu_factories(self) -> None:
        self.menu_factories = {
            "File": self._create_file_menu,
            "Edit": self._create_edit_menu,
            "Display": self._create_display_menu,
            "Mode": self._create_mode_menu,
            "Constraints": self._create_constraints_menu,
            "View": self._create_view_menu,
            "Help": self._create_help_menu,
        }

    def _build_native_menu_bar(self) -> None:
        self.native_menu_bar = tk.Menu(self.root, tearoff=False)
        self.root.configure(menu=self.native_menu_bar)
        self._refresh_native_menu_bar()

    def _refresh_native_menu_bar(self) -> None:
        if not IS_WINDOWS or not self.menu_factories:
            return
        self.native_menu_vars = []
        menu_bar = tk.Menu(self.root, tearoff=False)
        for name in ("File", "Edit", "Display", "Mode", "Constraints", "View", "Help"):
            submenu = tk.Menu(menu_bar, tearoff=False)
            self._populate_native_menu(submenu, self.menu_factories[name]())
            menu_bar.add_cascade(label=name, menu=submenu)
        self.native_menu_bar = menu_bar
        self.root.configure(menu=menu_bar)

    def _populate_native_menu(self, menu: tk.Menu, entries: list[tuple]) -> None:
        for entry in entries:
            kind = entry[0]
            if kind == "separator":
                menu.add_separator()
            elif kind == "command":
                menu.add_command(label=entry[1], command=entry[2], accelerator=entry[3] or "")
            elif kind == "disabled":
                menu.add_command(label=entry[1], state=tk.DISABLED, accelerator=entry[2] if len(entry) > 2 else "")
            elif kind == "submenu":
                submenu = tk.Menu(menu, tearoff=False)
                self._populate_native_menu(submenu, entry[2])
                menu.add_cascade(label=entry[1], menu=submenu)
            elif kind == "toggle":
                var = tk.BooleanVar(value=bool(entry[2]))
                self.native_menu_vars.append(var)
                menu.add_checkbutton(label=entry[1], variable=var, command=entry[3])

    def _build_menu_strip(self) -> None:
        self.menu_row = tk.Frame(self.main, bg="#e5e5e5", bd=1, relief=tk.RAISED)
        self.menu_row.pack(fill=tk.X)

        menu_names = ("File", "Edit", "Display", "Mode", "Constraints", "View", "Help")
        for name in menu_names:
            label = tk.Label(
                self.menu_row,
                text=name,
                bg="#e5e5e5",
                fg="#111827",
                padx=7,
                pady=2,
                bd=0,
            )
            label.pack(side=tk.LEFT)
            label.bind("<Enter>", lambda _event, menu_name=name: self._show_menu(menu_name))
            label.bind("<ButtonPress-1>", lambda _event, menu_name=name: self._toggle_menu(menu_name))
            self.menu_labels[name] = label

    def _build_toolbar(self) -> None:
        toolbar_bg = "#f3f4f6" if IS_WINDOWS else "#e5e5e5"
        self.toolbar_row = tk.Frame(
            self.main,
            bg=toolbar_bg,
            bd=0 if IS_WINDOWS else 1,
            relief=tk.FLAT if IS_WINDOWS else tk.RAISED,
        )
        self.toolbar_row.pack(fill=tk.X)
        self.toolbar_hint = tk.Label(
            self.main,
            bg="#ffffe1",
            fg="#111111",
            bd=1,
            relief=tk.SOLID,
            padx=5,
            pady=2,
            font=("TkDefaultFont", 9),
        )
        self.toolbar_hint.place_forget()
        self._add_toolbar_button("new_file", self._new_project, "normal")
        self._add_toolbar_button("open_file", self.open_file, "normal")
        self._add_toolbar_button("save", self.save_project, "disabled")
        self._add_toolbar_separator()
        self._add_toolbar_button("cut", self._cut_selection, "disabled")
        self._add_toolbar_button("copy", self._copy_selection, "disabled")
        self._add_toolbar_button("paste", self._paste_selection, "disabled")
        self._add_toolbar_button("delete", self._delete_selection, "disabled")
        self._add_toolbar_separator()
        self._add_toolbar_button("undo", self.undo, "disabled")
        self._add_toolbar_button("redo", self.redo, "disabled")
        self._add_toolbar_separator()
        self._add_toolbar_button("camera", self.render_project, "disabled")
        self._add_toolbar_separator()
        self._add_toolbar_button("grid", self._toggle_grid_visibility, "normal")
        self._add_toolbar_button("strand_number", self._toggle_strand_labels, "normal")
        self._add_toolbar_button("show_sticky_end", self._toggle_sticky_end_visibility, "normal")
        self._add_toolbar_separator()
        for mode, _label, icon_name, _shortcut in SELECTION_MODES:
            self._add_toolbar_button(icon_name, lambda selected_mode=mode: self._set_selection_mode(selected_mode), "normal")
        self._add_toolbar_separator()
        self._add_toolbar_button("create_helix", lambda: self._set_selection_mode("create_helix"), "normal")
        self._add_toolbar_button("create_base", lambda: self._set_selection_mode("create_base"), "normal")
        self._add_toolbar_separator()
        self._add_toolbar_button("move", lambda: self._set_selection_mode("move"), "disabled")
        self._add_toolbar_button("move_fixed", self._move_selection_fixed, "disabled")
        self._add_toolbar_button("rotate", lambda: self._set_selection_mode("rotate"), "disabled")
        self._add_toolbar_button("rotate_fixed", self._rotate_selection_fixed, "disabled")
        self._add_toolbar_separator()
        self._add_toolbar_button("set_sequence", self.set_sequence, "disabled")
        self._add_toolbar_button("reset_sequence", self.reset_sequence, "disabled")
        self._add_toolbar_button("generate_sequence", self.generate_sequences, "disabled")
        self._add_toolbar_separator()
        self._add_toolbar_button("change_backbone_color", self.change_backbone_color, "disabled")
        self._add_toolbar_button("reset_backbone_color", self.reset_backbone_color, "disabled")
        self._add_toolbar_separator()
        self._add_toolbar_button("create_complementary", self._create_complementary, "disabled")
        self._add_toolbar_button("delete_complementary", self._delete_complementary, "disabled")
        self._add_toolbar_button("create_down", self._create_down, "disabled")
        self._add_toolbar_button("delete_down", self._delete_down, "disabled")
        self._add_toolbar_button("create_sticky_end", self._create_sticky_end, "disabled")
        self._add_toolbar_button("delete_sticky_end", self._delete_sticky_end, "disabled")

    def _add_toolbar_button(self, icon_name: str, command, initial_state: str) -> None:
        image = self._load_toolbar_icon(icon_name)
        toolbar_bg = "#f3f4f6" if IS_WINDOWS else "#e5e5e5"
        slot = tk.Frame(
            self.toolbar_row,
            bg=toolbar_bg,
            width=32,
            height=32,
            highlightthickness=1,
            highlightbackground=toolbar_bg,
            highlightcolor=toolbar_bg,
            bd=0,
        )
        slot.pack(side=tk.LEFT, padx=1, pady=1)
        slot.pack_propagate(False)
        button = tk.Button(
            slot,
            image=image,
            command=command,
            state=initial_state,
            width=26,
            height=26,
            bg=toolbar_bg,
            activebackground=toolbar_bg if IS_WINDOWS else "#d4d4d4",
            relief=tk.FLAT,
            bd=0 if IS_WINDOWS else 1,
            highlightthickness=0,
            padx=1,
            pady=1,
        )
        button.pack(padx=2, pady=2)
        hint_text = TOOLBAR_HINTS.get(icon_name, icon_name.replace("_", " ").title())
        for widget in (slot, button):
            widget.bind("<Enter>", lambda _event, name=icon_name, text=hint_text: self._handle_toolbar_enter(name, text))
            widget.bind("<Leave>", lambda _event, name=icon_name: self._handle_toolbar_leave(name))
            widget.bind("<ButtonPress-1>", lambda _event, name=icon_name: self._handle_toolbar_press(name))
            widget.bind("<ButtonRelease-1>", lambda _event, name=icon_name: self._handle_toolbar_release(name))
        self.toolbar_buttons[icon_name] = button
        self.toolbar_slots[icon_name] = slot
        self._set_toolbar_button_enabled(icon_name, initial_state != "disabled")

    def _add_toolbar_separator(self) -> None:
        separator = tk.Frame(self.toolbar_row, bg="#cbd5e1" if IS_WINDOWS else "#808080", width=1, height=26)
        separator.pack(side=tk.LEFT, padx=4, pady=2)

    def _load_toolbar_icon(self, icon_name: str) -> tk.PhotoImage:
        if icon_name in self.toolbar_icons:
            return self.toolbar_icons[icon_name]
        icon_path = _resource_root() / "Toolbar_icons" / f"{icon_name}.png"
        if icon_path.exists():
            image = tk.PhotoImage(file=str(icon_path))
            max_dimension = max(image.width(), image.height(), 1)
            step = max(1, round(max_dimension / 18))
            scaled = image.subsample(step, step) if step > 1 else image
            self.toolbar_icons[icon_name] = scaled
            return scaled
        if icon_name == "select_base":
            image = self._build_select_base_icon()
            self.toolbar_icons[icon_name] = image
            return image
        if icon_name == "select_box":
            image = self._build_select_box_icon()
            self.toolbar_icons[icon_name] = image
            return image
        image = tk.PhotoImage(file=str(icon_path))
        max_dimension = max(image.width(), image.height(), 1)
        step = max(1, round(max_dimension / 18))
        scaled = image.subsample(step, step) if step > 1 else image
        self.toolbar_icons[icon_name] = scaled
        return scaled

    def _load_disabled_toolbar_icon(self, icon_name: str) -> tk.PhotoImage:
        if icon_name in self.toolbar_disabled_icons:
            return self.toolbar_disabled_icons[icon_name]
        disabled = _disabled_toolbar_icon(self._load_toolbar_icon(icon_name))
        self.toolbar_disabled_icons[icon_name] = disabled
        return disabled

    def _set_toolbar_button_enabled(self, icon_name: str, enabled: bool) -> None:
        button = self.toolbar_buttons.get(icon_name)
        if button is None:
            return
        button.configure(
            state=tk.NORMAL if enabled else tk.DISABLED,
            image=self._load_toolbar_icon(icon_name) if enabled else self._load_disabled_toolbar_icon(icon_name),
        )

    def _build_select_base_icon(self) -> tk.PhotoImage:
        image = tk.PhotoImage(width=18, height=18)
        for x in range(7, 11):
            for y in range(7, 11):
                image.put("#111111", (x, y))
        image.put("#111111", to=(8, 6, 10, 7))
        image.put("#111111", to=(8, 11, 10, 12))
        image.put("#111111", to=(6, 8, 7, 10))
        image.put("#111111", to=(11, 8, 12, 10))
        return image

    def _build_select_box_icon(self) -> tk.PhotoImage:
        image = tk.PhotoImage(width=18, height=18)
        left = 3
        top = 4
        right = 15
        bottom = 14
        for x in range(left, right + 1, 2):
            image.put("#111111", (x, top))
            image.put("#111111", (x, bottom))
        for y in range(top, bottom + 1, 2):
            image.put("#111111", (left, y))
            image.put("#111111", (right, y))
        return image

    def _create_file_menu(self) -> list[tuple]:
        entries: list[tuple] = [
            ("command", "New", self._new_project, _shortcut_label("new")),
            ("command", "Open...", self.open_file, _shortcut_label("open")),
        ]
        if self.project:
            entries.extend(
                [
                    ("command", "Save", self.save_project, _shortcut_label("save")),
                    ("command", "Save As...", self.save_json, _shortcut_label("save_as")),
                    ("separator",),
                    ("command", "Render...", self.render_project, None),
                    ("command", "Export Sequences...", self.export_sequences, None),
                ]
            )
        else:
            entries.extend(
                [
                    ("disabled", "Save"),
                    ("disabled", "Save As..."),
                    ("separator",),
                    ("disabled", "Render..."),
                    ("disabled", "Export Sequences..."),
                ]
            )
        recent_entries = self._create_recent_file_menu_entries()
        if recent_entries:
            entries.append(("separator",))
            entries.extend(recent_entries)
        entries.extend(
            [
                ("separator",),
                ("command", "Exit", self.root.destroy, None),
            ]
        )
        return entries

    def _create_edit_menu(self) -> list[tuple]:
        entries: list[tuple] = []
        if self.undo_stack:
            entries.append(("command", "Undo", self.undo, _shortcut_label("undo")))
        else:
            entries.append(("disabled", "Undo", _shortcut_label("undo")))
        if self.redo_stack:
            entries.append(("command", "Redo", self.redo, _shortcut_label("redo")))
        else:
            entries.append(("disabled", "Redo", _shortcut_label("redo")))
        entries.extend(
            [
                ("separator",),
                ("command", "Cut", self._cut_selection, _shortcut_label("cut")) if self._can_cut_selection() else ("disabled", "Cut", _shortcut_label("cut")),
                ("command", "Copy", self._copy_selection, _shortcut_label("copy")) if self._can_copy_selection() else ("disabled", "Copy", _shortcut_label("copy")),
                ("command", "Paste", self._paste_selection, _shortcut_label("paste")) if self._can_paste_selection() else ("disabled", "Paste", _shortcut_label("paste")),
            ]
        )
        if self.project and self.selected_indices:
            entries.append(("command", "Delete", self._delete_selection, "Del"))
        else:
            entries.append(("disabled", "Delete", "Del"))
        if self._can_create_complementary():
            entries.append(("command", "Create Complementary", self._create_complementary, None))
        else:
            entries.append(("disabled", "Create Complementary"))
        if self._can_delete_complementary():
            entries.append(("command", "Delete Complementary", self._delete_complementary, None))
        else:
            entries.append(("disabled", "Delete Complementary"))
        if self._can_create_down():
            entries.append(("command", "Create Down", self._create_down, None))
        else:
            entries.append(("disabled", "Create Down"))
        if self._can_delete_down():
            entries.append(("command", "Delete Down", self._delete_down, None))
        else:
            entries.append(("disabled", "Delete Down"))
        if self._can_create_sticky_end():
            entries.append(("command", "Create Sticky End", self._create_sticky_end, None))
        else:
            entries.append(("disabled", "Create Sticky End"))
        if self._can_delete_sticky_end():
            entries.append(("command", "Delete Sticky End", self._delete_sticky_end, None))
        else:
            entries.append(("disabled", "Delete Sticky End"))
        entries.append(("separator",))
        if self.project:
            entries.append(("command", "Select All", self._select_all, _shortcut_label("select_all")))
        else:
            entries.append(("disabled", "Select All", _shortcut_label("select_all")))
        entries.extend(
            [
                ("separator",),
                ("command", "Set Sequence", self.set_sequence, None) if self._can_set_sequence() else ("disabled", "Set Sequence"),
                ("command", "Reset Sequence", self.reset_sequence, None) if self._can_reset_sequence() else ("disabled", "Reset Sequence"),
                ("command", "Fill Sequence", self.fill_sequence, None) if self._can_fill_sequence() else ("disabled", "Fill Sequence"),
                ("command", "Generate Sequence", self.generate_sequences, None) if self.project else ("disabled", "Generate Sequence"),
            ]
        )
        return entries

    def _create_mode_menu(self) -> list[tuple]:
        entries = [
            ("command", label, lambda selected_mode=mode: self._set_selection_mode(selected_mode), shortcut)
            for mode, label, _icon_name, shortcut in SELECTION_MODES
        ]
        entries.append(("separator",))
        entries.extend(
            [
                ("command", label, lambda selected_mode=mode: self._set_selection_mode(selected_mode), shortcut)
                for mode, label, _icon_name, shortcut in CREATE_MODES
            ]
        )
        entries.append(("separator",))
        if self.project and self.selected_indices:
            entries.extend(
                [
                    ("command", "Move", lambda: self._set_selection_mode("move"), None),
                    ("command", "Rotate", lambda: self._set_selection_mode("rotate"), None),
                ]
            )
        else:
            entries.extend(
                [
                    ("disabled", "Move"),
                    ("disabled", "Rotate"),
                ]
            )
        return entries

    def _create_view_menu(self) -> list[tuple]:
        return [
            ("command", "Reset View Split", self._reset_split_layout, None),
            ("separator",),
            ("toggle", "Grid", self.show_grid, self._toggle_grid_visibility),
            ("toggle", "Strand Number", self.show_strand_labels, self._toggle_strand_labels),
            ("toggle", "Sticky Ends", self.show_sticky_end, self._toggle_sticky_end_visibility),
        ]

    def _create_display_menu(self) -> list[tuple]:
        return [
            (
                "submenu",
                "Nucleotide Size",
                self._create_size_choice_menu(
                    values=(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 8.0),
                    current=self.display_sphere_size,
                    setter=lambda value: self._set_display_sizes(sphere=value),
                ),
            ),
            (
                "submenu",
                "Backbone Size",
                self._create_size_choice_menu(
                    values=(1, 2, 3, 4, 5),
                    current=self.display_backbone_width,
                    setter=lambda value: self._set_display_sizes(backbone=int(value)),
                ),
            ),
            (
                "submenu",
                "Base Pair Size",
                self._create_size_choice_menu(
                    values=(1, 2, 3, 4, 5),
                    current=self.display_basepair_width,
                    setter=lambda value: self._set_display_sizes(basepair=int(value)),
                ),
            ),
            ("separator",),
            ("command", "Change Backbone Color...", self.change_backbone_color, None)
            if self._can_change_backbone_color()
            else ("disabled", "Change Backbone Color..."),
            ("command", "Reset Backbone Color", self.reset_backbone_color, None)
            if self._can_reset_backbone_color()
            else ("disabled", "Reset Backbone Color"),
            ("separator",),
            ("command", "Grid Extents...", self._edit_grid_extents, None),
            (
                "submenu",
                "Color Scheme",
                self._create_color_scheme_menu(),
            ),
        ]

    def _create_constraints_menu(self) -> list[tuple]:
        return [
            ("command", "View/Edit Constraints...", self._edit_constraints, None),
            ("command", "Customize Constraints...", self._customize_constraints, None),
        ]

    def _create_help_menu(self) -> list[tuple]:
        return [
            ("command", "About TiamatPy", self._show_about_dialog, None),
        ]

    def _create_placeholder_menu(self, name: str) -> list[tuple]:
        return [
            ("disabled", f"{name} actions will go here"),
        ]

    def _create_size_choice_menu(self, values: tuple[float | int, ...], current: float | int, setter) -> list[tuple]:
        entries: list[tuple] = []
        for value in values:
            numeric_value = float(value)
            label = f"{int(value) if numeric_value.is_integer() else numeric_value:g} Pixel"
            if numeric_value != 1.0:
                label += "s"
            marker = "* " if abs(numeric_value - float(current)) < 1e-6 else "  "
            entries.append(("command", f"{marker}{label}", lambda selected=value: setter(selected), None))
        return entries

    def _create_color_scheme_menu(self) -> list[tuple]:
        return [
            (
                "command",
                f"{'* ' if self.color_scheme == 'dark' else '  '}Dark Mode",
                lambda: self._set_color_scheme("dark"),
                None,
            ),
            (
                "command",
                f"{'* ' if self.color_scheme == 'bright' else '  '}Bright Mode",
                lambda: self._set_color_scheme("bright"),
                None,
            ),
        ]

    def _create_recent_file_menu_entries(self) -> list[tuple]:
        entries: list[tuple] = []
        for index, path in enumerate(self.recent_files[:4], start=1):
            entries.append(
                ("command", f"{index} {_display_recent_file_path(path)}", lambda selected_path=path: self._open_recent_file(selected_path), None)
            )
        return entries

    def _show_menu(self, name: str) -> None:
        if name not in self.menu_factories or name not in self.menu_labels:
            return
        if self._active_menu_name == name:
            return
        self._close_active_menu()
        label = self.menu_labels[name]
        self._set_menu_label_state(name, active=True)
        self._active_menu_name = name
        x = label.winfo_x()
        y = self.menu_row.winfo_y() + label.winfo_y() + label.winfo_height() - 2
        self.menu_popup = self._build_menu_popup(self.menu_factories[name](), is_submenu=False)
        self.menu_popup.place(x=x, y=y)

    def _toggle_menu(self, name: str) -> None:
        if self._active_menu_name == name:
            self._close_active_menu()
            return
        self._show_menu(name)

    def _close_active_menu(self) -> None:
        if self._active_menu_name is None:
            return
        name = self._active_menu_name
        self._active_menu_name = None
        self._close_active_submenu()
        if self.menu_popup is not None:
            self.menu_popup.destroy()
            self.menu_popup = None
        self._set_menu_label_state(name, active=False)

    def _set_menu_label_state(self, name: str, active: bool) -> None:
        label = self.menu_labels.get(name)
        if label is None:
            return
        if active:
            label.configure(bg="#2563eb", fg="#ffffff")
        else:
            label.configure(bg="#e5e5e5", fg="#111827")

    def _close_menu_on_click(self, event: tk.Event) -> None:
        if self._active_menu_name is None:
            return
        widget = event.widget
        if self.menu_popup is not None and _widget_is_descendant(widget, self.menu_popup):
            return
        if self.submenu_popup is not None and _widget_is_descendant(widget, self.submenu_popup):
            return
        if widget in self.menu_labels.values():
            return
        self._close_active_menu()

    def _close_menu_on_escape(self, _event: tk.Event) -> str | None:
        if self._active_menu_name is None:
            return None
        self._close_active_menu()
        return "break"

    def _reset_split_layout(self) -> None:
        self._col_ratio = 0.5
        self._row_ratio = 0.5
        self._layout_viewports()

    def _set_display_sizes(
        self,
        sphere: float | None = None,
        backbone: int | None = None,
        basepair: int | None = None,
    ) -> None:
        if sphere is not None:
            self.display_sphere_size = sphere
        if backbone is not None:
            self.display_backbone_width = backbone
        if basepair is not None:
            self.display_basepair_width = basepair
        self._refresh_all()

    def _show_about_dialog(self) -> None:
        messagebox.showinfo(
            "About TiamatPy",
            "Tiamat Python version is programmed by Liangxiao Chen with help of Codex from OpenAI.\n\n"
            "Copyright (c) 2026 Liangxiao Chen. All rights reserved.",
        )

    def _edit_constraints(self) -> None:
        values = self._prompt_constraints_dialog("View/Edit Constraints", self.constraints_values)
        if values is None:
            return
        self._store_constraints_values(values)

    def _customize_constraints(self) -> None:
        values = self._prompt_constraints_dialog(
            "Customize Constraints",
            self.custom_constraints_values,
            CUSTOM_CONSTRAINT_COLUMNS,
        )
        if values is None:
            return
        self._store_custom_constraints_values(values)

    def _store_constraints_values(self, values: dict[str, dict[str, str]]) -> None:
        normalized = _normalize_constraints_values(values, fallback=values)
        if self.project is None:
            self.constraints_values = normalized
            return
        snapshot = self._capture_history_snapshot()
        self.constraints_values = normalized
        self.project.metadata["constraints"] = _copy_constraints_values(normalized)
        self._commit_history_snapshot(snapshot)
        self._refresh_all()

    def _store_custom_constraints_values(self, values: dict[str, dict[str, str]]) -> None:
        normalized = _normalize_constraints_values(values, columns=CUSTOM_CONSTRAINT_COLUMNS, fallback=values)
        if self.project is None:
            self.custom_constraints_values = normalized
            return
        snapshot = self._capture_history_snapshot()
        self.custom_constraints_values = normalized
        self.project.metadata["custom_constraints"] = _copy_constraints_values(normalized)
        self._commit_history_snapshot(snapshot)
        self._refresh_all()

    def _prompt_constraints_dialog(
        self,
        title: str,
        initial_values: dict[str, dict[str, str]],
        columns: tuple[tuple[str, str], ...] = CONSTRAINT_COLUMNS,
    ) -> dict[str, dict[str, str]] | None:
        dialog = tk.Toplevel(self.root)
        dialog.withdraw()
        dialog.title(title)
        dialog.configure(bg="#e5e5e5")
        dialog.resizable(False, False)
        dialog.transient(self.root)

        values = {
            row_key: {
                column_key: tk.StringVar(value=initial_values[row_key][column_key])
                for column_key, _column_label in columns
            }
            for row_key, _row_label in CONSTRAINT_ROWS
        }
        result: dict[str, dict[str, str]] | None = None
        reset_values = _default_constraints_values() if columns == CONSTRAINT_COLUMNS else None

        self._bind_helicity_constraint_traces(values, columns)

        body = tk.Frame(dialog, bg="#e5e5e5", padx=12, pady=12)
        body.pack(fill=tk.BOTH, expand=True)
        for column_index in range(1, len(columns) + 1):
            body.grid_columnconfigure(column_index, minsize=230, uniform="constraint_columns")

        for column_index, (_column_key, column_label) in enumerate(columns, start=1):
            tk.Label(
                body,
                text=column_label,
                bg="#e5e5e5",
                anchor="center",
                justify=tk.CENTER,
            ).grid(row=0, column=column_index, padx=6, pady=(0, 10))

        first_entry = None
        for row_index, (row_key, row_label) in enumerate(CONSTRAINT_ROWS, start=1):
            tk.Label(body, text=row_label, bg="#e5e5e5", anchor="w").grid(row=row_index, column=0, sticky="w", padx=(0, 10), pady=4)
            for column_index, (column_key, _column_label) in enumerate(columns, start=1):
                entry = tk.Entry(body, textvariable=values[row_key][column_key], width=12)
                entry.grid(row=row_index, column=column_index, sticky="ew", padx=6, pady=4)
                if first_entry is None:
                    first_entry = entry

        buttons = tk.Frame(body, bg="#e5e5e5")
        buttons.grid(row=len(CONSTRAINT_ROWS) + 1, column=1, columnspan=max(len(columns), 1), sticky="e", pady=(12, 0))

        def reset_defaults() -> None:
            if reset_values is None:
                return
            for row_key, _row_label in CONSTRAINT_ROWS:
                for column_key, _column_label in columns:
                    values[row_key][column_key].set(reset_values[row_key][column_key])

        def confirm() -> None:
            nonlocal result
            try:
                parsed = {
                    row_key: {
                        column_key: values[row_key][column_key].get().strip()
                        for column_key, _column_label in columns
                    }
                    for row_key, _row_label in CONSTRAINT_ROWS
                }
                result = _validate_constraints_values(parsed, columns=columns)
            except ValueError as exc:
                messagebox.showerror(title, str(exc), parent=dialog)
                return
            dialog.destroy()

        column = 0
        if reset_values is not None:
            reset_button = tk.Button(buttons, text="Reset", width=9, command=reset_defaults, highlightthickness=0)
            reset_button.grid(row=0, column=column, sticky="ew", padx=(0, 10))
            column += 1
        ok_button = tk.Button(buttons, text="OK", width=9, command=confirm, highlightthickness=0)
        ok_button.grid(row=0, column=column, sticky="ew", padx=(0, 10))
        column += 1
        cancel_button = tk.Button(buttons, text="Cancel", width=9, command=dialog.destroy, highlightthickness=0)
        cancel_button.grid(row=0, column=column, sticky="ew")

        dialog.bind("<Return>", lambda _event: confirm())
        dialog.bind("<Escape>", lambda _event: dialog.destroy())
        dialog.protocol("WM_DELETE_WINDOW", dialog.destroy)

        self.root.update_idletasks()
        dialog.update_idletasks()
        x = self.root.winfo_rootx() + max((self.root.winfo_width() - dialog.winfo_width()) // 2, 0)
        y = self.root.winfo_rooty() + max((self.root.winfo_height() - dialog.winfo_height()) // 2, 0)
        dialog.geometry(f"+{x}+{y}")
        dialog.deiconify()
        dialog.lift()
        dialog.grab_set()
        if first_entry is not None:
            first_entry.focus_set()
            first_entry.selection_range(0, tk.END)
        dialog.wait_window()
        return result

    def _bind_helicity_constraint_traces(
        self,
        values: dict[str, dict[str, tk.StringVar]],
        columns: tuple[tuple[str, str], ...],
    ) -> None:
        if "rotation_per_base_pair_deg" not in values or "helicity_bp_per_turn" not in values:
            return
        for column_key, _column_label in columns:
            rotation_var = values["rotation_per_base_pair_deg"][column_key]
            helicity_var = values["helicity_bp_per_turn"][column_key]
            state = {"updating": False}

            def sync_from_rotation(*_args: object, column_key: str = column_key) -> None:
                if state["updating"]:
                    return
                text = values["rotation_per_base_pair_deg"][column_key].get().strip()
                if text == "":
                    return
                try:
                    derived = _format_helicity_value(_helicity_from_rotation_value(float(text)))
                except ValueError:
                    return
                if values["helicity_bp_per_turn"][column_key].get().strip() == derived:
                    return
                state["updating"] = True
                values["helicity_bp_per_turn"][column_key].set(derived)
                state["updating"] = False

            def sync_from_helicity(*_args: object, column_key: str = column_key) -> None:
                if state["updating"]:
                    return
                text = values["helicity_bp_per_turn"][column_key].get().strip()
                if text == "":
                    return
                try:
                    derived = _format_rotation_value(_rotation_from_helicity_value(float(text)))
                except ValueError:
                    return
                if values["rotation_per_base_pair_deg"][column_key].get().strip() == derived:
                    return
                state["updating"] = True
                values["rotation_per_base_pair_deg"][column_key].set(derived)
                state["updating"] = False

            rotation_var.trace_add("write", sync_from_rotation)
            helicity_var.trace_add("write", sync_from_helicity)

    def _prompt_create_helix_dialog(self, initial_count: int) -> dict[str, str | int | float] | None:
        dialog = tk.Toplevel(self.root)
        dialog.withdraw()
        dialog.title("Create Helix")
        dialog.configure(bg="#e5e5e5")
        dialog.resizable(False, False)
        dialog.transient(self.root)

        number_var = tk.StringVar(value=str(max(1, initial_count)))
        rotation_var = tk.StringVar(value="0")
        molecule_var = tk.StringVar(value=CREATE_MOLECULE_OPTIONS[0])
        family_options = list(CREATE_HELIX_FAMILY_OPTIONS)
        if _custom_constraints_defined(self.custom_constraints_values):
            family_options.append("X-type")
        structure_var = tk.StringVar(value="B-type")
        type_var = tk.StringVar(value="Helix")
        result: dict[str, str | int | float] | None = None

        body = tk.Frame(dialog, bg="#e5e5e5", padx=12, pady=12)
        body.pack(fill=tk.BOTH, expand=True)

        tk.Label(body, text="Number of Bases:", bg="#e5e5e5", anchor="w").grid(row=0, column=0, sticky="w", pady=4, padx=(0, 8))
        first_entry = tk.Entry(body, textvariable=number_var, width=18)
        first_entry.grid(row=0, column=1, sticky="ew", pady=4)

        tk.Label(body, text="Backbone Rotation:", bg="#e5e5e5", anchor="w").grid(row=1, column=0, sticky="w", pady=4, padx=(0, 8))
        rotation_entry = tk.Entry(body, textvariable=rotation_var, width=18)
        rotation_entry.grid(row=1, column=1, sticky="ew", pady=4)

        tk.Label(body, text="Molecule:", bg="#e5e5e5", anchor="w").grid(row=2, column=0, sticky="nw", pady=(8, 4), padx=(0, 8))
        molecule_frame = tk.Frame(body, bg="#e5e5e5")
        molecule_frame.grid(row=2, column=1, sticky="w", pady=(8, 4))
        for row, label in enumerate(CREATE_MOLECULE_OPTIONS):
            radio = tk.Radiobutton(
                molecule_frame,
                text=label,
                value=label,
                variable=molecule_var,
                bg="#e5e5e5",
                activebackground="#e5e5e5",
                highlightthickness=0,
                anchor="w",
                justify=tk.LEFT,
                selectcolor="#ffffff",
            )
            radio.grid(row=row // 2, column=row % 2, sticky="w", padx=(0, 18), pady=2)

        tk.Label(body, text="Structure:", bg="#e5e5e5", anchor="w").grid(row=3, column=0, sticky="nw", pady=(8, 4), padx=(0, 8))
        structure_frame = tk.Frame(body, bg="#e5e5e5")
        structure_frame.grid(row=3, column=1, sticky="w", pady=(8, 4))
        for row, label in enumerate(family_options):
            radio = tk.Radiobutton(
                structure_frame,
                text=label,
                value=label,
                variable=structure_var,
                bg="#e5e5e5",
                activebackground="#e5e5e5",
                highlightthickness=0,
                anchor="w",
                justify=tk.LEFT,
                selectcolor="#ffffff",
            )
            radio.grid(row=row // 2, column=row % 2, sticky="w", padx=(0, 18), pady=2)

        tk.Label(body, text="Type:", bg="#e5e5e5", anchor="w").grid(row=4, column=0, sticky="nw", pady=(8, 4), padx=(0, 8))
        type_frame = tk.Frame(body, bg="#e5e5e5")
        type_frame.grid(row=4, column=1, sticky="w", pady=(8, 4))
        initial_focus = None
        for row, label in enumerate(CREATE_STRUCTURE_TYPES):
            radio = tk.Radiobutton(
                type_frame,
                text=label,
                value=label,
                variable=type_var,
                bg="#e5e5e5",
                activebackground="#e5e5e5",
                highlightthickness=0,
                anchor="w",
                justify=tk.LEFT,
                selectcolor="#ffffff",
            )
            radio.grid(row=row // 2, column=row % 2, sticky="w", padx=(0, 18), pady=2)
            if initial_focus is None:
                initial_focus = radio

        preview = self._load_toolbar_icon("create_helix")
        preview_label = tk.Label(body, image=preview, bg="#e5e5e5")
        preview_label.image = preview
        preview_label.grid(row=0, column=2, rowspan=5, sticky="n", padx=(20, 0))

        buttons = tk.Frame(body, bg="#e5e5e5")
        buttons.grid(row=0, column=3, rowspan=5, sticky="ne", padx=(14, 0))

        def confirm() -> None:
            nonlocal result
            try:
                count = int(number_var.get().strip())
            except ValueError:
                messagebox.showerror("Create Helix", "Number of Bases must be an integer.", parent=dialog)
                return
            if count <= 0:
                messagebox.showerror("Create Helix", "Number of Bases must be greater than 0.", parent=dialog)
                return
            try:
                backbone_rotation = float(rotation_var.get().strip() or "0")
            except ValueError:
                messagebox.showerror("Create Helix", "Backbone Rotation must be numeric.", parent=dialog)
                return
            result = {
                "count": count,
                "backbone_rotation": backbone_rotation,
                "molecule": molecule_var.get(),
                "structure_family": structure_var.get(),
                "structure_type": type_var.get(),
            }
            dialog.destroy()

        ok_button = tk.Button(buttons, text="OK", width=9, command=confirm, highlightthickness=0)
        ok_button.grid(row=0, column=0, sticky="ew", pady=(0, 8))
        cancel_button = tk.Button(buttons, text="Cancel", width=9, command=dialog.destroy, highlightthickness=0)
        cancel_button.grid(row=1, column=0, sticky="ew")

        dialog.bind("<Return>", lambda _event: confirm())
        dialog.bind("<Escape>", lambda _event: dialog.destroy())
        dialog.protocol("WM_DELETE_WINDOW", dialog.destroy)

        self.root.update_idletasks()
        dialog.update_idletasks()
        x = self.root.winfo_rootx() + max((self.root.winfo_width() - dialog.winfo_width()) // 2, 0)
        y = self.root.winfo_rooty() + max((self.root.winfo_height() - dialog.winfo_height()) // 2, 0)
        dialog.geometry(f"+{x}+{y}")
        dialog.deiconify()
        dialog.lift()
        dialog.grab_set()
        first_entry.focus_set()
        first_entry.selection_range(0, tk.END)
        if initial_focus is not None:
            initial_focus.configure(takefocus=0)
        dialog.wait_window()
        return result

    def _create_helix_from_interaction(self, key: str, start_canvas: tuple[int, int], end_canvas: tuple[int, int], dragged: bool) -> None:
        state = self.viewports[key]
        width = max(1, state.canvas.winfo_width())
        height = max(1, state.canvas.winfo_height())
        center = state.scene_center or (_project_center(self.project) if self.project and self.project.bases else (0.0, 0.0, 0.0))
        scale = _view_scale(self.project, self.grid_extents, state, width, height, center)
        start_point = _canvas_to_world(
            start_canvas[0],
            start_canvas[1],
            state.camera,
            width,
            height,
            state.mode,
            center,
            scale,
        )
        if dragged:
            end_point = _canvas_to_world(
                end_canvas[0],
                end_canvas[1],
                state.camera,
                width,
                height,
                state.mode,
                center,
                scale,
            )
            direction = (
                end_point[0] - start_point[0],
                end_point[1] - start_point[1],
                end_point[2] - start_point[2],
            )
            axis_length = sqrt(direction[0] * direction[0] + direction[1] * direction[1] + direction[2] * direction[2])
        else:
            direction = (1.0, 0.0, 0.0)
            axis_length = 0.0

        estimated_count = _estimate_created_base_count(axis_length, self.constraints_values)
        options = self._prompt_create_helix_dialog(estimated_count)
        if options is None:
            return
        structure_family = str(options["structure_family"])
        active_constraints = self.custom_constraints_values if structure_family == "X-type" else self.constraints_values

        snapshot = self._capture_history_snapshot()
        if self.project is None:
            self.project = TiamatProject(
                bases={},
                metadata={
                    "constraints": _copy_constraints_values(self.constraints_values),
                    "custom_constraints": _copy_constraints_values(self.custom_constraints_values),
                },
            )
        try:
            created_indices = _append_created_structure(
                self.project,
                start_point=start_point,
                direction=direction,
                count=int(options["count"]),
                backbone_rotation=float(options["backbone_rotation"]),
                molecule_option=_create_helix_molecule_option(str(options["molecule"]), structure_family),
                structure_type=str(options["structure_type"]),
                constraints_values=active_constraints,
            )
        except ValueError as exc:
            messagebox.showerror("Create Helix", str(exc), parent=self.root)
            if snapshot == self._capture_history_snapshot():
                return
            self._restore_history_snapshot(snapshot)
            self._refresh_all()
            return

        self.project.metadata["constraints"] = _copy_constraints_values(self.constraints_values)
        self.project.metadata["custom_constraints"] = _copy_constraints_values(self.custom_constraints_values)
        self.selected_indices = set(created_indices)
        self.primary_selected = created_indices[0] if created_indices else None
        for viewport_state in self.viewports.values():
            _mark_viewport_for_fit(viewport_state, recenter=True)
        self._commit_history_snapshot(snapshot)
        self._refresh_all()

    def _create_helix_drag_endpoint(self, key: str, state: ViewportState, event: tk.Event) -> tuple[int, int]:
        current = (event.x, event.y)
        start = state.selection_box_start
        if start is None or not _has_shift_modifier(event.state):
            return current
        width = max(1, state.canvas.winfo_width())
        height = max(1, state.canvas.winfo_height())
        center = state.scene_center or (_project_center(self.project) if self.project and self.project.bases else (0.0, 0.0, 0.0))
        scale = _view_scale(self.project, self.grid_extents, state, width, height, center)
        return _snap_create_helix_endpoint(
            start,
            current,
            state.camera,
            width,
            height,
            state.mode,
            center,
            scale,
        )

    def _prompt_free_strand_length_dialog(self, default_count: int) -> int | None:
        dialog = tk.Toplevel(self.root)
        dialog.withdraw()
        dialog.title("Create Free Strand")
        dialog.configure(bg="#e5e5e5")
        dialog.resizable(False, False)
        dialog.transient(self.root)

        count_var = tk.StringVar(value=str(max(1, default_count)))
        result: dict[str, int | None] = {"value": None}

        body = tk.Frame(dialog, bg="#e5e5e5", padx=12, pady=12)
        body.pack(fill=tk.BOTH, expand=True)

        tk.Label(body, text="Free strand length:", bg="#e5e5e5", anchor="w").grid(row=0, column=0, sticky="w", padx=(0, 8))
        entry = tk.Entry(body, textvariable=count_var, width=12)
        entry.grid(row=0, column=1, sticky="ew")

        buttons = tk.Frame(body, bg="#e5e5e5")
        buttons.grid(row=0, column=2, padx=(12, 0))

        def confirm() -> None:
            try:
                count = int(count_var.get().strip())
            except ValueError:
                messagebox.showerror("Create Free Strand", "Free strand length must be an integer.", parent=dialog)
                return
            if count <= 0:
                messagebox.showerror("Create Free Strand", "Free strand length must be greater than 0.", parent=dialog)
                return
            result["value"] = count
            dialog.destroy()

        tk.Button(buttons, text="OK", width=9, command=confirm, highlightthickness=0).grid(row=0, column=0, sticky="ew", pady=(0, 8))
        tk.Button(buttons, text="Cancel", width=9, command=dialog.destroy, highlightthickness=0).grid(row=1, column=0, sticky="ew")

        dialog.bind("<Return>", lambda _event: confirm())
        dialog.bind("<Escape>", lambda _event: dialog.destroy())
        dialog.protocol("WM_DELETE_WINDOW", dialog.destroy)

        self.root.update_idletasks()
        dialog.update_idletasks()
        x = self.root.winfo_rootx() + max((self.root.winfo_width() - dialog.winfo_width()) // 2, 0)
        y = self.root.winfo_rooty() + max((self.root.winfo_height() - dialog.winfo_height()) // 2, 0)
        dialog.geometry(f"+{x}+{y}")
        dialog.deiconify()
        dialog.lift()
        dialog.grab_set()
        entry.focus_set()
        entry.selection_range(0, tk.END)
        dialog.wait_window()
        return result["value"]

    def _free_strand_click_target(self, key: str, x: int, y: int) -> tuple[FreeStrandEndpoint | None, str | None]:
        state = self.viewports[key]
        best_index = _nearest_projected_index(state.projected, x, y)
        if self.project and best_index is not None:
            base = self.project.bases[best_index]
            roles = _base_terminal_roles(base)
            if not roles:
                return None, "Free strand endpoints must be empty space or terminal 3'/5' bases."
            return FreeStrandEndpoint(point=base.position, base_index=best_index, terminal_roles=roles, molecule=base.molecule), None
        width = max(1, state.canvas.winfo_width())
        height = max(1, state.canvas.winfo_height())
        center = state.scene_center or (_project_center(self.project) if self.project and self.project.bases else (0.0, 0.0, 0.0))
        scale = _view_scale(self.project, self.grid_extents, state, width, height, center)
        point = _canvas_to_world(x, y, state.camera, width, height, state.mode, center, scale)
        return FreeStrandEndpoint(point=point), None

    def _handle_create_free_strand_click(self, key: str, event: tk.Event) -> None:
        endpoint, error = self._free_strand_click_target(key, event.x, event.y)
        if endpoint is None:
            if error:
                messagebox.showerror("Create Free Strand", error, parent=self.root)
            return
        if self.pending_free_strand_start is None:
            self.pending_free_strand_start = endpoint
            self.pending_free_strand_hover_view = key
            self.pending_free_strand_hover_canvas = (event.x, event.y)
            if endpoint.base_index is not None:
                self.selected_indices = {endpoint.base_index}
                self.primary_selected = endpoint.base_index
            else:
                self.selected_indices.clear()
                self.primary_selected = None
            self._refresh_all()
            return

        try:
            plan = _resolve_free_strand_creation(self.pending_free_strand_start, endpoint)
        except ValueError as exc:
            messagebox.showerror("Create Free Strand", str(exc), parent=self.root)
            return

        default_count = _estimate_free_strand_base_count(self.pending_free_strand_start.point, endpoint.point)
        count = self._prompt_free_strand_length_dialog(default_count)
        if count is None:
            return

        snapshot = self._capture_history_snapshot()
        if self.project is None:
            self.project = TiamatProject(
                bases={},
                metadata={
                    "constraints": _copy_constraints_values(self.constraints_values),
                    "custom_constraints": _copy_constraints_values(self.custom_constraints_values),
                },
            )
        try:
            created_indices = _append_free_strand(
                self.project,
                start_point=plan["start_point"],
                end_point=plan["end_point"],
                count=count,
                molecule=plan["molecule"],
                start_anchor=plan["start_anchor"],
                end_anchor=plan["end_anchor"],
            )
        except ValueError as exc:
            messagebox.showerror("Create Free Strand", str(exc), parent=self.root)
            if snapshot != self._capture_history_snapshot():
                self._restore_history_snapshot(snapshot)
                self._refresh_all()
            return

        self.project.metadata["constraints"] = _copy_constraints_values(self.constraints_values)
        self.project.metadata["custom_constraints"] = _copy_constraints_values(self.custom_constraints_values)
        self._clear_pending_free_strand()
        self.selected_indices = set(created_indices)
        self.primary_selected = created_indices[0] if created_indices else None
        for viewport_state in self.viewports.values():
            _mark_viewport_for_fit(viewport_state, recenter=True)
        self._commit_history_snapshot(snapshot)
        self._refresh_all()

    def _move_selection_fixed(self) -> None:
        if not self.project or not self.selected_indices:
            return
        values = self._prompt_vector_dialog("Move Selected Bases", "nm")
        if values is None:
            return
        self._apply_translation_to_selection(*values)

    def _rotate_selection_fixed(self) -> None:
        if not self.project or not self.selected_indices:
            return
        values = self._prompt_vector_dialog("Rotate Selected Bases", "deg")
        if values is None:
            return
        rotation = _matrix_multiply(
            _rotation_matrix_z(values[2]),
            _matrix_multiply(_rotation_matrix_y(values[1]), _rotation_matrix_x(values[0])),
        )
        self._apply_rotation_to_selection(rotation)

    def _prompt_vector_dialog(self, title: str, unit_label: str) -> tuple[float, float, float] | None:
        dialog = tk.Toplevel(self.root)
        dialog.title(title)
        dialog.resizable(False, False)
        dialog.transient(self.root)
        dialog.grab_set()

        values = {
            "X": tk.StringVar(value="0"),
            "Y": tk.StringVar(value="0"),
            "Z": tk.StringVar(value="0"),
        }
        result: dict[str, tuple[float, float, float] | None] = {"value": None}

        body = ttk.Frame(dialog, padding=12)
        body.grid(row=0, column=0, sticky="nsew")

        for row, axis in enumerate(("X", "Y", "Z")):
            ttk.Label(body, text=f"{axis}:").grid(row=row, column=0, sticky="w", padx=(0, 6), pady=4)
            entry = ttk.Entry(body, textvariable=values[axis], width=10)
            entry.grid(row=row, column=1, sticky="ew", pady=4)
            ttk.Label(body, text=unit_label).grid(row=row, column=2, sticky="w", padx=(6, 0), pady=4)
            if row == 0:
                entry.focus_set()
                entry.selection_range(0, tk.END)

        buttons = ttk.Frame(body)
        buttons.grid(row=0, column=3, rowspan=3, sticky="ns", padx=(12, 0))
        buttons.columnconfigure(0, weight=1)
        ok_button = ttk.Button(buttons, text="OK")
        ok_button.grid(row=0, column=0, sticky="ew", pady=(0, 6))
        cancel_button = ttk.Button(buttons, text="Cancel", command=dialog.destroy)
        cancel_button.grid(row=1, column=0, sticky="ew")

        def confirm() -> None:
            try:
                parsed = tuple(float(values[axis].get().strip() or "0") for axis in ("X", "Y", "Z"))
            except ValueError:
                messagebox.showerror("Invalid value", f"Enter numeric {unit_label} values.", parent=dialog)
                return
            result["value"] = parsed
            dialog.destroy()

        ok_button.configure(command=confirm)
        dialog.bind("<Return>", lambda _event: confirm())
        dialog.bind("<Escape>", lambda _event: dialog.destroy())
        dialog.wait_window()
        return result["value"]

    def _apply_translation_to_selection(self, dx: float, dy: float, dz: float) -> None:
        if not self.project or not self.selected_indices:
            return
        snapshot = self._capture_history_snapshot()
        for index in self.selected_indices:
            base = self.project.bases[index]
            base.x += dx
            base.y += dy
            base.z += dz
        self._commit_history_snapshot(snapshot)
        self._refresh_all()

    def _apply_rotation_to_selection(
        self,
        rotation: tuple[tuple[float, float, float], tuple[float, float, float], tuple[float, float, float]],
    ) -> None:
        if not self.project or not self.selected_indices:
            return
        snapshot = self._capture_history_snapshot()
        center = _selection_center(self.project, self.selected_indices)
        original_positions = {
            index: self.project.bases[index].position
            for index in self.selected_indices
        }
        _apply_rotation_to_project(self.project, original_positions, center, rotation)
        self._commit_history_snapshot(snapshot)
        self._refresh_all()

    def _show_toolbar_hint(self, widget: tk.Widget, text: str) -> None:
        if self.toolbar_hint is None:
            return
        self.toolbar_hint.configure(text=text)
        self.main.update_idletasks()
        hint_width = self.toolbar_hint.winfo_reqwidth()
        hint_height = self.toolbar_hint.winfo_reqheight()
        x = widget.winfo_rootx() - self.main.winfo_rootx() + widget.winfo_width() - 6
        y = widget.winfo_rooty() - self.main.winfo_rooty() + widget.winfo_height() - 3
        max_x = max(0, self.main.winfo_width() - hint_width - 4)
        max_y = max(0, self.main.winfo_height() - hint_height - 4)
        self.toolbar_hint.place(x=min(x, max_x), y=min(y, max_y))
        self.toolbar_hint.lift()

    def _hide_toolbar_hint(self) -> None:
        if self.toolbar_hint is not None:
            self.toolbar_hint.place_forget()

    def _set_selection_mode(self, mode: str) -> None:
        if mode not in INTERACTION_MODE_LABELS:
            return
        if mode in TRANSFORM_MODE_NAMES and (not self.project or not self.selected_indices):
            return
        if self.selection_mode == "create_base":
            self._clear_pending_free_strand()
        self.selection_mode = mode
        for state in self.viewports.values():
            state.selection_box_start = None
            state.selection_box_end = None
            state.drag_positions = None
            state.drag_center = None
            state.drag_projected_center = None
        self._hide_toolbar_hint()
        self._refresh_all()

    def _handle_toolbar_enter(self, icon_name: str, hint_text: str) -> None:
        button = self.toolbar_buttons.get(icon_name)
        if button is None:
            return
        self._set_toolbar_button_visual(icon_name, hover=True)
        self._show_toolbar_hint(button, hint_text)

    def _handle_toolbar_leave(self, icon_name: str) -> None:
        self._set_toolbar_button_visual(icon_name, hover=False)
        self._hide_toolbar_hint()

    def _handle_toolbar_press(self, icon_name: str) -> None:
        self._hide_toolbar_hint()
        self._set_toolbar_button_visual(icon_name, hover=True, pressed=True)

    def _handle_toolbar_release(self, icon_name: str) -> None:
        self._set_toolbar_button_visual(icon_name, hover=True)

    def _set_toolbar_button_visual(self, icon_name: str, hover: bool, pressed: bool = False) -> None:
        button = self.toolbar_buttons.get(icon_name)
        slot = self.toolbar_slots.get(icon_name)
        if button is None or slot is None:
            return
        active = INTERACTION_MODE_ICONS.get(self.selection_mode) == icon_name
        if icon_name == "strand_number":
            active = self.show_strand_labels
        if icon_name == "show_sticky_end":
            active = self.show_sticky_end
        if icon_name == "grid":
            active = self.show_grid
        base_bg = "#f3f4f6" if IS_WINDOWS else "#e5e5e5"
        hover_outline = "#60a5fa" if IS_WINDOWS else "#2563eb"
        active_bg = "#dbeafe" if IS_WINDOWS else "#dcdcdc"
        pressed_bg = "#bfdbfe" if IS_WINDOWS else "#dcdcdc"
        outline = hover_outline if hover or active else base_bg
        if IS_WINDOWS:
            slot_bg = pressed_bg if pressed and str(button.cget("state")) != tk.DISABLED else active_bg if active else base_bg
            slot.configure(
                bg=slot_bg,
                highlightbackground=outline,
                highlightcolor=outline,
            )
            button.configure(
                relief=tk.FLAT,
                bg=slot_bg,
                activebackground=slot_bg,
                bd=0,
            )
            return
        slot.configure(bg=base_bg, highlightbackground=outline, highlightcolor=outline)
        if pressed and str(button.cget("state")) != tk.DISABLED:
            button.configure(relief=tk.SUNKEN, bg=pressed_bg)
            button.pack_configure(padx=(3, 1), pady=(3, 1))
            return
        if active and str(button.cget("state")) != tk.DISABLED:
            button.configure(relief=tk.SUNKEN, bg=active_bg)
            button.pack_configure(padx=(3, 1), pady=(3, 1))
            return
        button.configure(relief=tk.FLAT, bg=base_bg)
        button.pack_configure(padx=2, pady=2)

    def _toggle_strand_labels(self) -> None:
        self.show_strand_labels = not self.show_strand_labels
        self._refresh_all()

    def _toggle_sticky_end_visibility(self) -> None:
        self.show_sticky_end = not self.show_sticky_end
        self._refresh_all()

    def _toggle_grid_visibility(self) -> None:
        self.show_grid = not self.show_grid
        self._refresh_all()

    def _theme(self) -> dict[str, str]:
        return COLOR_SCHEMES.get(self.color_scheme, COLOR_SCHEMES["dark"])

    def _apply_color_scheme(self) -> None:
        theme = self._theme()
        self.viewport_host.configure(background=theme["viewport_bg"])
        self.vertical_sash.configure(background=theme["sash_bg"])
        self.horizontal_sash.configure(background=theme["sash_bg"])
        for state in self.viewports.values():
            state.canvas.master.configure(background=theme["viewport_bg"])
            state.canvas.configure(
                background=theme["canvas_bg"],
                highlightbackground=theme["canvas_border"],
            )

    def _set_color_scheme(self, scheme: str) -> None:
        if scheme not in COLOR_SCHEMES or self.color_scheme == scheme:
            return
        self.color_scheme = scheme
        self._apply_color_scheme()
        self._refresh_all()

    def _edit_grid_extents(self) -> None:
        dialog = tk.Toplevel(self.root)
        dialog.withdraw()
        dialog.title("Grid Extents")
        dialog.configure(bg="#e5e5e5")
        dialog.resizable(False, False)
        dialog.transient(self.root)

        labels = (
            ("X min (nm)", "x_min"),
            ("X max (nm)", "x_max"),
            ("Y min (nm)", "y_min"),
            ("Y max (nm)", "y_max"),
            ("Z min (nm)", "z_min"),
            ("Z max (nm)", "z_max"),
            ("Spacing (nm)", "spacing"),
        )
        values = {key: tk.StringVar(value=f"{self.grid_extents[key]:g}") for _label, key in labels}
        result = {"accepted": False}

        body = tk.Frame(dialog, bg="#e5e5e5", padx=12, pady=12)
        body.pack(fill=tk.BOTH, expand=True)

        for row, (label_text, key) in enumerate(labels):
            tk.Label(body, text=label_text, bg="#e5e5e5", anchor="w").grid(row=row, column=0, sticky="w", pady=3)
            tk.Entry(body, textvariable=values[key], width=12).grid(row=row, column=1, sticky="ew", padx=(10, 0), pady=3)

        buttons = tk.Frame(body, bg="#e5e5e5")
        buttons.grid(row=0, column=2, rowspan=len(labels), sticky="ne", padx=(14, 0))

        def confirm() -> None:
            try:
                parsed = {key: float(values[key].get().strip()) for _label, key in labels}
            except ValueError:
                messagebox.showerror("Invalid value", "Enter numeric grid extents and spacing.", parent=dialog)
                return
            if parsed["spacing"] <= 0.0:
                messagebox.showerror("Invalid value", "Grid spacing must be greater than 0.", parent=dialog)
                return
            for axis in ("x", "y", "z"):
                if parsed[f"{axis}_min"] >= parsed[f"{axis}_max"]:
                    messagebox.showerror(
                        "Invalid value",
                        f"{axis.upper()} min must be smaller than {axis.upper()} max.",
                        parent=dialog,
                    )
                    return
            self.grid_extents = parsed
            result["accepted"] = True
            dialog.destroy()

        ok_button = tk.Button(buttons, text="OK", width=9, command=confirm, highlightthickness=0)
        ok_button.grid(row=0, column=0, sticky="ew", pady=(0, 8))
        cancel_button = tk.Button(buttons, text="Cancel", width=9, command=dialog.destroy, highlightthickness=0)
        cancel_button.grid(row=1, column=0, sticky="ew")

        dialog.bind("<Return>", lambda _event: confirm())
        dialog.bind("<Escape>", lambda _event: dialog.destroy())
        dialog.protocol("WM_DELETE_WINDOW", dialog.destroy)

        self.root.update_idletasks()
        dialog.update_idletasks()
        x = self.root.winfo_rootx() + max((self.root.winfo_width() - dialog.winfo_width()) // 2, 0)
        y = self.root.winfo_rooty() + max((self.root.winfo_height() - dialog.winfo_height()) // 2, 0)
        dialog.geometry(f"+{x}+{y}")
        dialog.deiconify()
        dialog.lift()
        dialog.grab_set()
        dialog.wait_window()

        if result["accepted"]:
            self._refresh_all()

    def _capture_history_snapshot(self) -> dict[str, object]:
        return _snapshot_history_state(self.project, self.selected_indices, self.primary_selected)

    def _can_create_complementary(self) -> bool:
        return bool(self.project and self.project.can_create_complementary(self.selected_indices))

    def _can_delete_complementary(self) -> bool:
        return bool(self.project and self.project.can_delete_complementary(self.selected_indices))

    def _can_create_down(self) -> bool:
        return bool(self.project and self.project.can_create_down(self.selected_indices))

    def _can_delete_down(self) -> bool:
        return bool(self.project and self.project.can_delete_down(self.selected_indices))

    def _can_create_sticky_end(self) -> bool:
        return bool(self.project and self.project.can_create_sticky_end(self.selected_indices))

    def _can_delete_sticky_end(self) -> bool:
        return bool(self.project and self.project.can_delete_sticky_end(self.selected_indices))

    def _selected_strand(self):
        if not self.project:
            return None
        return _selected_strand_for_indices(self.project, self.selected_indices)

    def _can_set_sequence(self) -> bool:
        return self._selected_strand() is not None

    def _can_reset_sequence(self) -> bool:
        return bool(self.project and self.selected_indices)

    def _can_copy_selection(self) -> bool:
        return bool(self.project and self.selected_indices)

    def _can_cut_selection(self) -> bool:
        return self._can_copy_selection()

    def _can_paste_selection(self) -> bool:
        return self.structure_clipboard is not None and bool(self.structure_clipboard.get("bases"))

    def _can_change_backbone_color(self) -> bool:
        return bool(self.project and self.selected_indices)

    def _can_reset_backbone_color(self) -> bool:
        return bool(self.project and self.selected_indices)

    def _can_fill_sequence(self) -> bool:
        return bool(
            self.project
            and any(
                base.nucleotide is None and (base.across is None or base.across not in self.project.bases)
                for base in self.project.bases.values()
            )
        )

    def _create_complementary(self) -> None:
        if not self.project:
            return
        snapshot = self._capture_history_snapshot()
        if self.project.create_complementary(self.selected_indices):
            self._commit_history_snapshot(snapshot)
            self._refresh_all()

    def _copy_selection(self) -> None:
        if not self._can_copy_selection() or not self.project:
            return
        self.structure_clipboard = _selection_clipboard_payload(self.project, self.selected_indices)
        self.structure_clipboard_paste_count = 0
        self._update_toolbar_state()
        self._update_status()

    def _cut_selection(self) -> None:
        if not self._can_cut_selection():
            return
        self._copy_selection()
        self._delete_selection()

    def _paste_selection(self) -> None:
        if not self._can_paste_selection():
            return
        snapshot = self._capture_history_snapshot()
        self.structure_clipboard_paste_count += 1
        offset = _clipboard_paste_offset(self.structure_clipboard_paste_count)
        self.project, new_indices = _paste_structure_clipboard(self.project, self.structure_clipboard, offset)
        self.selected_indices = set(new_indices)
        self.primary_selected = min(new_indices) if new_indices else None
        self._commit_history_snapshot(snapshot)
        self._refresh_all()

    def _delete_complementary(self) -> None:
        if not self.project:
            return
        snapshot = self._capture_history_snapshot()
        if self.project.delete_complementary(self.selected_indices):
            self._commit_history_snapshot(snapshot)
            self._refresh_all()

    def _create_down(self) -> None:
        if not self.project:
            return
        snapshot = self._capture_history_snapshot()
        if self.project.create_down(self.selected_indices):
            self._commit_history_snapshot(snapshot)
            self._refresh_all()

    def _delete_down(self) -> None:
        if not self.project:
            return
        snapshot = self._capture_history_snapshot()
        if self.project.delete_down(self.selected_indices):
            self._commit_history_snapshot(snapshot)
            self._refresh_all()

    def _create_sticky_end(self) -> None:
        if not self.project:
            return
        snapshot = self._capture_history_snapshot()
        created, error = self.project.create_sticky_end(self.selected_indices)
        if not created:
            if error:
                messagebox.showerror("Create Sticky End", error, parent=self.root)
            return
        self._commit_history_snapshot(snapshot)
        self._refresh_all()

    def _delete_sticky_end(self) -> None:
        if not self.project:
            return
        snapshot = self._capture_history_snapshot()
        if self.project.delete_sticky_end(self.selected_indices):
            self._commit_history_snapshot(snapshot)
            self._refresh_all()

    def change_backbone_color(self) -> None:
        if not self.project or not self.selected_indices:
            return
        initial_color = _initial_backbone_color(self.project, self.selected_indices)
        _rgb, selected_color = colorchooser.askcolor(
            color=initial_color,
            parent=self.root,
            title="Change Backbone Color",
        )
        if not selected_color:
            return
        snapshot = self._capture_history_snapshot()
        if self.project.set_backbone_color(self.selected_indices, selected_color):
            self._commit_history_snapshot(snapshot)
            self._refresh_all()

    def reset_backbone_color(self) -> None:
        if not self.project or not self.selected_indices:
            return
        snapshot = self._capture_history_snapshot()
        if self.project.reset_backbone_color(self.selected_indices):
            self._commit_history_snapshot(snapshot)
            self._refresh_all()

    def set_sequence(self) -> None:
        strand = self._selected_strand()
        if not self.project or strand is None:
            return

        first_base = self.project.bases[strand.base_indices[0]]
        molecule = str(first_base.molecule).upper()
        legend_items = _sequence_legend_items(molecule)
        allowed_bases = _allowed_sequence_bases(molecule)
        current_sequence = strand.sequence(self.project)

        dialog = tk.Toplevel(self.root)
        dialog.withdraw()
        dialog.title("Set Sequence")
        dialog.configure(bg="#e5e5e5")
        dialog.resizable(False, False)
        dialog.transient(self.root)

        body = tk.Frame(dialog, bg="#e5e5e5", padx=12, pady=12)
        body.pack(fill=tk.BOTH, expand=True)

        legend = tk.Frame(body, bg="#e5e5e5")
        legend.grid(row=0, column=0, sticky="nw", padx=(0, 16))
        for row, text in enumerate(legend_items):
            tk.Label(
                legend,
                text=text,
                bg="#e5e5e5",
                anchor="w",
                justify=tk.LEFT,
            ).grid(row=row, column=0, sticky="w")

        editor_frame = tk.Frame(body, bg="#e5e5e5")
        editor_frame.grid(row=0, column=1, sticky="nsew")
        body.grid_columnconfigure(1, weight=1)

        sequence_text = tk.Text(
            editor_frame,
            width=38,
            height=3,
            wrap="word",
            highlightthickness=1,
            highlightbackground="#9ca3af",
            highlightcolor="#2563eb",
            bd=0,
        )
        sequence_text.grid(row=0, column=0, sticky="ew")
        sequence_text.insert("1.0", current_sequence)

        helper_var = tk.StringVar(value=f"Must be exactly {len(strand.base_indices)} nucleotides long.")
        current_var = tk.StringVar(value=f"Current length: {len(current_sequence)} nucleotides.")
        error_var = tk.StringVar(value="")

        tk.Label(
            editor_frame,
            textvariable=helper_var,
            bg="#e5e5e5",
            anchor="w",
            justify=tk.LEFT,
        ).grid(row=1, column=0, sticky="w", pady=(8, 0))
        tk.Label(
            editor_frame,
            textvariable=current_var,
            bg="#e5e5e5",
            anchor="w",
            justify=tk.LEFT,
        ).grid(row=2, column=0, sticky="w", pady=(4, 0))
        error_label = tk.Label(
            editor_frame,
            textvariable=error_var,
            bg="#e5e5e5",
            fg="#b91c1c",
            anchor="w",
            justify=tk.LEFT,
        )
        error_label.grid(row=3, column=0, sticky="w", pady=(4, 0))

        buttons = tk.Frame(body, bg="#e5e5e5")
        buttons.grid(row=0, column=2, sticky="ne", padx=(14, 0))
        result: dict[str, str | None] = {"sequence": None}

        ok_button = tk.Button(buttons, text="OK", width=9, highlightthickness=0)
        ok_button.grid(row=0, column=0, sticky="ew", pady=(0, 8))
        cancel_button = tk.Button(buttons, text="Cancel", width=9, command=dialog.destroy, highlightthickness=0)
        cancel_button.grid(row=1, column=0, sticky="ew")

        def update_state() -> tuple[str, list[str], bool]:
            raw_text = sequence_text.get("1.0", "end-1c")
            normalized, invalid = _normalize_sequence_input(raw_text, allowed_bases)
            is_valid = not invalid and len(normalized) == len(strand.base_indices)
            current_var.set(f"Current length: {len(normalized)} nucleotides.")
            if invalid:
                error_var.set(f"Invalid characters: {' '.join(invalid)}")
            elif len(normalized) != len(strand.base_indices):
                error_var.set("")
            else:
                error_var.set("")
            ok_button.configure(state=tk.NORMAL if is_valid else tk.DISABLED)
            return normalized, invalid, is_valid

        def schedule_update(_event=None) -> None:
            dialog.after_idle(update_state)

        def confirm() -> None:
            normalized, _invalid, is_valid = update_state()
            if not is_valid:
                return
            result["sequence"] = normalized
            dialog.destroy()

        ok_button.configure(command=confirm)
        sequence_text.bind("<KeyRelease>", schedule_update)
        sequence_text.bind("<<Paste>>", schedule_update)
        sequence_text.bind("<<Cut>>", schedule_update)
        dialog.bind("<Return>", lambda _event: confirm())
        dialog.bind("<Escape>", lambda _event: dialog.destroy())
        dialog.protocol("WM_DELETE_WINDOW", dialog.destroy)

        update_state()
        self.root.update_idletasks()
        dialog.update_idletasks()
        x = self.root.winfo_rootx() + max((self.root.winfo_width() - dialog.winfo_width()) // 2, 0)
        y = self.root.winfo_rooty() + max((self.root.winfo_height() - dialog.winfo_height()) // 2, 0)
        dialog.geometry(f"+{x}+{y}")
        dialog.deiconify()
        dialog.lift()
        dialog.grab_set()
        sequence_text.focus_set()
        sequence_text.tag_add("sel", "1.0", "end-1c")
        dialog.wait_window()

        normalized_sequence = result["sequence"]
        if not normalized_sequence:
            return

        snapshot = self._capture_history_snapshot()
        for index, nucleotide in zip(strand.base_indices, normalized_sequence):
            self.project.set_nucleotide(index, None if nucleotide == "N" else nucleotide)
        self._commit_history_snapshot(snapshot)
        self._refresh_all()

    def reset_sequence(self) -> None:
        if not self.project or not self.selected_indices:
            return
        snapshot = self._capture_history_snapshot()
        for index in sorted(self.selected_indices):
            if index in self.project.bases:
                self.project.set_nucleotide(index, None)
        self._commit_history_snapshot(snapshot)
        self._refresh_all()

    def fill_sequence(self) -> None:
        if not self.project or not self._can_fill_sequence():
            return

        dialog = tk.Toplevel(self.root)
        dialog.withdraw()
        dialog.title("Fill Sequence")
        dialog.configure(bg="#e5e5e5")
        dialog.resizable(False, False)
        dialog.transient(self.root)

        value_var = tk.StringVar(value="A")
        result: dict[str, str | None] = {"value": None}

        body = tk.Frame(dialog, bg="#e5e5e5", padx=12, pady=12)
        body.pack(fill=tk.BOTH, expand=True)

        tk.Label(
            body,
            text="Fill the unpaired generic base to",
            bg="#e5e5e5",
            anchor="w",
        ).grid(row=0, column=0, sticky="w", padx=(0, 18))

        options_frame = tk.Frame(body, bg="#e5e5e5")
        options_frame.grid(row=1, column=0, sticky="nw", pady=(8, 0))

        initial_focus = None
        for row, (label, value) in enumerate((("A", "A"), ("T/U", "T"), ("C", "C"), ("G", "G"))):
            radio = tk.Radiobutton(
                options_frame,
                text=label,
                value=value,
                variable=value_var,
                bg="#e5e5e5",
                activebackground="#e5e5e5",
                highlightthickness=0,
                anchor="w",
                justify=tk.LEFT,
                selectcolor="#ffffff",
            )
            radio.grid(row=row, column=0, sticky="w", pady=2)
            if initial_focus is None:
                initial_focus = radio

        buttons = tk.Frame(body, bg="#e5e5e5")
        buttons.grid(row=1, column=1, sticky="ne", padx=(14, 0))

        def confirm() -> None:
            result["value"] = value_var.get()
            dialog.destroy()

        ok_button = tk.Button(buttons, text="OK", width=9, command=confirm, highlightthickness=0)
        ok_button.grid(row=0, column=0, sticky="ew", pady=(0, 8))
        cancel_button = tk.Button(buttons, text="Cancel", width=9, command=dialog.destroy, highlightthickness=0)
        cancel_button.grid(row=1, column=0, sticky="ew")

        dialog.bind("<Return>", lambda _event: confirm())
        dialog.bind("<Escape>", lambda _event: dialog.destroy())
        dialog.protocol("WM_DELETE_WINDOW", dialog.destroy)

        self.root.update_idletasks()
        dialog.update_idletasks()
        x = self.root.winfo_rootx() + max((self.root.winfo_width() - dialog.winfo_width()) // 2, 0)
        y = self.root.winfo_rooty() + max((self.root.winfo_height() - dialog.winfo_height()) // 2, 0)
        dialog.geometry(f"+{x}+{y}")
        dialog.deiconify()
        dialog.lift()
        dialog.grab_set()
        if initial_focus is not None:
            initial_focus.focus_set()
        dialog.wait_window()

        fill_value = result["value"]
        if fill_value is None:
            return

        snapshot = self._capture_history_snapshot()
        changed = self.project.fill_unpaired_generic_bases(fill_value)
        if changed:
            self._commit_history_snapshot(snapshot)
            self._refresh_all()

    def _commit_history_snapshot(self, snapshot: dict[str, object] | None) -> None:
        if snapshot is None:
            return
        if snapshot == self._capture_history_snapshot():
            return
        self.undo_stack.append(snapshot)
        if len(self.undo_stack) > MAX_HISTORY_STEPS:
            self.undo_stack = self.undo_stack[-MAX_HISTORY_STEPS:]
        self.redo_stack.clear()

    def _clear_history(self) -> None:
        self.undo_stack.clear()
        self.redo_stack.clear()

    def _restore_history_snapshot(self, snapshot: dict[str, object]) -> None:
        self.project, self.selected_indices, self.primary_selected = _restore_history_state(snapshot)
        if not self.project:
            self.constraints_values = _default_constraints_values()
            self.custom_constraints_values = _blank_custom_constraints_values()
            self._reset_viewports()
            return
        self.constraints_values = _project_constraints(self.project)
        self.custom_constraints_values = _project_custom_constraints(self.project)
        valid_indices = set(self.project.bases)
        self.selected_indices.intersection_update(valid_indices)
        if self.primary_selected not in self.selected_indices:
            self.primary_selected = min(self.selected_indices) if self.selected_indices else None

    def undo(self) -> None:
        if not self.undo_stack:
            return
        current = self._capture_history_snapshot()
        snapshot = self.undo_stack.pop()
        self.redo_stack.append(current)
        if len(self.redo_stack) > MAX_HISTORY_STEPS:
            self.redo_stack = self.redo_stack[-MAX_HISTORY_STEPS:]
        self._restore_history_snapshot(snapshot)
        self._refresh_all()

    def redo(self) -> None:
        if not self.redo_stack:
            return
        current = self._capture_history_snapshot()
        snapshot = self.redo_stack.pop()
        self.undo_stack.append(current)
        if len(self.undo_stack) > MAX_HISTORY_STEPS:
            self.undo_stack = self.undo_stack[-MAX_HISTORY_STEPS:]
        self._restore_history_snapshot(snapshot)
        self._refresh_all()

    def _load_recent_files(self) -> list[str]:
        try:
            payload = json.loads(VIEWER_STATE_PATH.read_text(encoding="utf-8"))
        except (FileNotFoundError, json.JSONDecodeError, OSError):
            return []
        paths = payload.get("recent_files")
        if not isinstance(paths, list):
            return []
        recent: list[str] = []
        for item in paths:
            path = str(item)
            if path and Path(path).suffix.lower() in OPENABLE_SUFFIXES and path not in recent:
                recent.append(path)
        return recent[:4]

    def _save_recent_files(self) -> None:
        try:
            VIEWER_STATE_PATH.write_text(
                json.dumps({"recent_files": self.recent_files[:4]}, indent=2) + "\n",
                encoding="utf-8",
            )
        except OSError:
            pass

    def _remember_recent_file(self, path: str | Path) -> None:
        normalized = str(Path(path).expanduser().resolve())
        if Path(normalized).suffix.lower() not in OPENABLE_SUFFIXES:
            return
        self.recent_files = [item for item in self.recent_files if item != normalized]
        self.recent_files.insert(0, normalized)
        self.recent_files = self.recent_files[:4]
        self._save_recent_files()

    def _open_project_path(self, path: str | Path) -> None:
        file_path = Path(path)
        self.project = load_project(file_path)
        self.constraints_values = _project_constraints(self.project)
        self.custom_constraints_values = _project_custom_constraints(self.project)
        self._remember_recent_file(file_path)
        self.selected_indices.clear()
        self.primary_selected = None
        self._clear_pending_free_strand()
        self._clear_history()
        self._reset_viewports()
        self._refresh_all()

    def _open_recent_file(self, path: str) -> None:
        try:
            self._open_project_path(path)
        except Exception as exc:
            self.recent_files = [item for item in self.recent_files if item != path]
            self._save_recent_files()
            messagebox.showerror("Open failed", str(exc))

    def _project_save_path(self) -> Path | None:
        if not self.project:
            return None
        source_path = self.project.metadata.get("source_path")
        if not source_path:
            return None
        path = Path(str(source_path))
        if path.suffix.lower() not in {".json", ".dnajson"}:
            return None
        return path

    def _write_project_json(self, path: str | Path) -> None:
        if not self.project:
            return
        self.project.metadata["constraints"] = _copy_constraints_values(self.constraints_values)
        self.project.metadata["custom_constraints"] = _copy_constraints_values(self.custom_constraints_values)
        save_json_project(self.project, path)
        self.project.metadata["source_path"] = str(Path(path))
        self.project.metadata["source_format"] = "dnajson"
        self._remember_recent_file(path)
        self._refresh_all()

    def _build_menu_popup(self, entries: list[tuple], is_submenu: bool) -> tk.Frame:
        popup = tk.Frame(
            self.main,
            bg="#dcdcdc",
            bd=1,
            relief=tk.SOLID,
            highlightthickness=1,
            highlightbackground="#7a7a7a",
        )
        max_label = max((len(entry[1]) for entry in entries if entry[0] in {"command", "disabled", "submenu", "toggle"}), default=10)
        max_accel = max(
            (
                len(entry[3] or "")
                if entry[0] == "command"
                else len(entry[2] or "")
                if entry[0] == "disabled" and len(entry) > 2
                else 1
                for entry in entries
                if entry[0] in {"command", "submenu", "disabled"}
            ),
            default=0,
        )

        row_index = 0
        for entry in entries:
            kind = entry[0]
            if kind == "separator":
                line = tk.Frame(popup, bg="#9a9a9a", height=1)
                line.grid(row=row_index, column=0, columnspan=3, sticky="ew", padx=2, pady=3)
                row_index += 1
                continue

            text = entry[1]
            accelerator = ""
            command = None
            submenu_entries = None
            mark_text = ""
            enabled = kind in {"command", "submenu", "toggle"}
            if kind == "command":
                command = entry[2]
                accelerator = entry[3] or ""
            elif kind == "disabled":
                accelerator = entry[2] if len(entry) > 2 else ""
            elif kind == "toggle":
                mark_text = _menu_toggle_mark(bool(entry[2]))
                command = entry[3]
            elif kind == "submenu":
                submenu_entries = entry[2]
                accelerator = ">"

            row = tk.Frame(popup, bg="#dcdcdc")
            row.grid(row=row_index, column=0, columnspan=3, sticky="ew")
            row.grid_columnconfigure(0, minsize=20)
            row.grid_columnconfigure(1, minsize=max_label * 8 + 8)
            row.grid_columnconfigure(2, minsize=max_accel * 7 + 22)

            mark = tk.Label(
                row,
                text=mark_text,
                bg="#dcdcdc",
                fg="#111111" if enabled else "#8f8f8f",
                anchor="w",
                padx=4,
                pady=3,
            )
            mark.grid(row=0, column=0, sticky="w")

            label = tk.Label(
                row,
                text=text,
                bg="#dcdcdc",
                fg="#111111" if enabled else "#8f8f8f",
                anchor="w",
                padx=4,
                pady=3,
            )
            label.grid(row=0, column=1, sticky="w")

            accel = tk.Label(
                row,
                text=accelerator,
                bg="#dcdcdc",
                fg="#111111" if enabled else "#8f8f8f",
                anchor="e",
                padx=8,
                pady=3,
            )
            accel.grid(row=0, column=2, sticky="e")

            widgets = (row, mark, label, accel)
            if kind == "submenu" and submenu_entries is not None:
                for widget in widgets:
                    widget.bind(
                        "<Enter>",
                        lambda _event, menu_row=row, row_widgets=widgets, submenu=submenu_entries: self._open_submenu(menu_row, row_widgets, submenu),
                    )
                    widget.bind(
                        "<Leave>",
                        lambda _event, row_widgets=widgets: self._restore_popup_row_state(row_widgets),
                    )
                    widget.bind(
                        "<ButtonRelease-1>",
                        lambda _event, menu_row=row, row_widgets=widgets, submenu=submenu_entries: self._open_submenu(menu_row, row_widgets, submenu),
                    )
            elif enabled and command is not None:
                for widget in widgets:
                    if is_submenu:
                        widget.bind("<Enter>", lambda _event, row_widgets=widgets: self._set_popup_row_state(row_widgets, True))
                    else:
                        widget.bind("<Enter>", lambda _event, row_widgets=widgets: self._activate_popup_row(row_widgets))
                    widget.bind("<Leave>", lambda _event, row_widgets=widgets: self._restore_popup_row_state(row_widgets))
                    widget.bind("<ButtonRelease-1>", lambda _event, callback=command: self._invoke_menu_command(callback))
            else:
                for widget in widgets:
                    widget.bind("<Enter>", lambda _event: self._close_active_submenu())

            row_index += 1

        return popup

    def _set_popup_row_state(self, widgets: tuple[tk.Widget, ...], active: bool) -> None:
        bg = "#0a64d8" if active else "#dcdcdc"
        fg = "#ffffff" if active else "#111111"
        for widget in widgets:
            if isinstance(widget, tk.Label):
                if widget.cget("fg") == "#8f8f8f":
                    continue
                widget.configure(bg=bg, fg=fg)
            else:
                widget.configure(bg=bg)

    def _activate_popup_row(self, widgets: tuple[tk.Widget, ...]) -> None:
        self._close_active_submenu()
        self._set_popup_row_state(widgets, True)

    def _restore_popup_row_state(self, widgets: tuple[tk.Widget, ...]) -> None:
        if self._active_submenu_widgets == widgets:
            return
        self._set_popup_row_state(widgets, False)

    def _open_submenu(self, row: tk.Frame, widgets: tuple[tk.Widget, ...], entries: list[tuple]) -> None:
        if self._active_submenu_widgets is not None and self._active_submenu_widgets != widgets:
            self._set_popup_row_state(self._active_submenu_widgets, False)
        self._close_active_submenu(reset_widgets=False)
        self._set_popup_row_state(widgets, True)
        self._active_submenu_widgets = widgets
        self.submenu_popup = self._build_menu_popup(entries, is_submenu=True)
        self.main.update_idletasks()
        x = row.winfo_rootx() - self.main.winfo_rootx() + row.winfo_width() - 1
        y = row.winfo_rooty() - self.main.winfo_rooty() - 1
        self.submenu_popup.place(x=x, y=y)

    def _close_active_submenu(self, reset_widgets: bool = True) -> None:
        if reset_widgets and self._active_submenu_widgets is not None:
            self._set_popup_row_state(self._active_submenu_widgets, False)
        self._active_submenu_widgets = None
        if self.submenu_popup is not None:
            self.submenu_popup.destroy()
            self.submenu_popup = None

    def _invoke_menu_command(self, command) -> None:
        self._hide_toolbar_hint()
        self._close_active_menu()
        command()

    def _bind_events(self) -> None:
        for key, state in self.viewports.items():
            state.canvas.bind("<Configure>", lambda _event, view_key=key: self._redraw_view(view_key))
            state.canvas.bind("<Motion>", lambda event, view_key=key: self._hover_view(view_key, event))
            state.canvas.bind("<Leave>", lambda _event, view_key=key: self._leave_view(view_key))
            state.canvas.bind("<ButtonPress-1>", lambda event, view_key=key: self._begin_drag(view_key, 1, event))
            state.canvas.bind("<B1-Motion>", lambda event, view_key=key: self._drag_view(view_key, 1, event))
            state.canvas.bind("<ButtonRelease-1>", lambda event, view_key=key: self._end_drag(view_key, 1, event))
            state.canvas.bind("<Double-Button-1>", lambda event, view_key=key: self._queue_base_editor_from_click(view_key, event))

            for button in (2, 3):
                state.canvas.bind(
                    f"<ButtonPress-{button}>",
                    lambda event, view_key=key, drag_button=button: self._begin_drag(view_key, drag_button, event),
                )
                state.canvas.bind(
                    f"<B{button}-Motion>",
                    lambda event, view_key=key, drag_button=button: self._drag_view(view_key, drag_button, event),
                )
                state.canvas.bind(
                    f"<ButtonRelease-{button}>",
                    lambda event, view_key=key, drag_button=button: self._end_drag(view_key, drag_button, event),
                )

            state.canvas.bind("<MouseWheel>", lambda event, view_key=key: self._on_wheel(view_key, event))
            state.canvas.bind("<Button-4>", lambda _event, view_key=key: self._on_scroll_zoom(view_key, 0.1))
            state.canvas.bind("<Button-5>", lambda _event, view_key=key: self._on_scroll_zoom(view_key, -0.1))

    def _hover_view(self, key: str, event: tk.Event) -> None:
        if self.selection_mode != "create_base" or self.pending_free_strand_start is None:
            return
        previous_view = self.pending_free_strand_hover_view
        self.pending_free_strand_hover_view = key
        self.pending_free_strand_hover_canvas = (event.x, event.y)
        if previous_view is not None and previous_view != key:
            self._redraw_view(previous_view)
        self._redraw_view(key)

    def _leave_view(self, key: str) -> None:
        if self.pending_free_strand_hover_view != key:
            return
        self.pending_free_strand_hover_view = None
        self.pending_free_strand_hover_canvas = None
        self._redraw_view(key)

    def _clear_pending_free_strand(self) -> None:
        self.pending_free_strand_start = None
        self.pending_free_strand_hover_view = None
        self.pending_free_strand_hover_canvas = None

        self.root.bind_all("<KeyPress>", self._on_keypress, add="+")
        self.root.bind_all("<Command-a>", self._select_all_event, add="+")
        self.root.bind_all("<Control-a>", self._select_all_event, add="+")
        self.root.bind_all("<Command-n>", self._new_project_event, add="+")
        self.root.bind_all("<Control-n>", self._new_project_event, add="+")
        self.root.bind_all("<Command-o>", self._open_file_event, add="+")
        self.root.bind_all("<Control-o>", self._open_file_event, add="+")
        self.root.bind_all("<Command-s>", self._save_project_event, add="+")
        self.root.bind_all("<Control-s>", self._save_project_event, add="+")
        self.root.bind_all("<Command-S>", self._save_project_as_event, add="+")
        self.root.bind_all("<Control-S>", self._save_project_as_event, add="+")
        self.root.bind_all("<Command-x>", self._cut_event, add="+")
        self.root.bind_all("<Control-x>", self._cut_event, add="+")
        self.root.bind_all("<Command-c>", self._copy_event, add="+")
        self.root.bind_all("<Control-c>", self._copy_event, add="+")
        self.root.bind_all("<Command-v>", self._paste_event, add="+")
        self.root.bind_all("<Control-v>", self._paste_event, add="+")
        self.root.bind_all("<Command-z>", self._undo_event, add="+")
        self.root.bind_all("<Control-z>", self._undo_event, add="+")
        self.root.bind_all("<Command-Z>", self._redo_event, add="+")
        self.root.bind_all("<Control-Z>", self._redo_event, add="+")
        self.root.bind_all("<ButtonPress-1>", self._close_menu_on_click, add="+")
        self.root.bind_all("<Escape>", self._close_menu_on_escape, add="+")

    def _layout_viewports(self) -> None:
        width = max(1, self.viewport_host.winfo_width())
        height = max(1, self.viewport_host.winfo_height())
        available_width = max(1, width - self._sash_thickness)
        available_height = max(1, height - self._sash_thickness)

        # Allow panes to collapse almost completely while keeping the sash visible.
        min_pane = 0
        left_width = int(available_width * self._col_ratio)
        top_height = int(available_height * self._row_ratio)
        left_width = max(min_pane, min(left_width, available_width - min_pane))
        top_height = max(min_pane, min(top_height, available_height - min_pane))
        right_width = available_width - left_width
        bottom_height = available_height - top_height

        positions = {
            "xz": (0, 0, left_width, top_height),
            "view3d": (left_width + self._sash_thickness, 0, right_width, top_height),
            "yz": (0, top_height + self._sash_thickness, left_width, bottom_height),
            "xy": (left_width + self._sash_thickness, top_height + self._sash_thickness, right_width, bottom_height),
        }
        for key, state in self.viewports.items():
            x, y, pane_width, pane_height = positions[key]
            state.canvas.master.place(x=x, y=y, width=pane_width, height=pane_height)

        self.vertical_sash.place(x=left_width, y=0, width=self._sash_thickness, height=height)
        self.horizontal_sash.place(x=0, y=top_height, width=width, height=self._sash_thickness)

    def _active_sash_set(self, name: str | None) -> None:
        self._active_sash = name

    def _drag_vertical_sash(self, event: tk.Event) -> None:
        if self._active_sash != "vertical":
            return
        width = max(1, self.viewport_host.winfo_width())
        available_width = max(1, width - self._sash_thickness)
        host_x = self.viewport_host.winfo_rootx()
        pointer_x = event.x_root - host_x
        ratio = (pointer_x - self._sash_thickness / 2.0) / available_width
        self._col_ratio = max(0.0, min(ratio, 1.0))
        self._layout_viewports()

    def _drag_horizontal_sash(self, event: tk.Event) -> None:
        if self._active_sash != "horizontal":
            return
        height = max(1, self.viewport_host.winfo_height())
        available_height = max(1, height - self._sash_thickness)
        host_y = self.viewport_host.winfo_rooty()
        pointer_y = event.y_root - host_y
        ratio = (pointer_y - self._sash_thickness / 2.0) / available_height
        self._row_ratio = max(0.0, min(ratio, 1.0))
        self._layout_viewports()

    def _refresh_all(self) -> None:
        self._update_toolbar_state()
        self._refresh_mode_toolbar_buttons()
        self._refresh_native_menu_bar()
        self._update_status()
        self._refresh_views_only()

    def _refresh_views_only(self) -> None:
        for key in self.viewports:
            self._redraw_view(key)

    def _refresh_mode_toolbar_buttons(self) -> None:
        for _mode, _label, icon_name, _shortcut in INTERACTION_MODES:
            self._set_toolbar_button_visual(icon_name, hover=False)
        self._set_toolbar_button_visual("strand_number", hover=False)
        self._set_toolbar_button_visual("show_sticky_end", hover=False)
        self._set_toolbar_button_visual("grid", hover=False)

    def _update_toolbar_state(self) -> None:
        has_project = self.project is not None
        has_selection = bool(self.selected_indices)
        if not has_selection and self.selection_mode in TRANSFORM_MODE_NAMES:
            self.selection_mode = "base"
        self._set_toolbar_button_enabled("save", has_project)
        self._set_toolbar_button_enabled("camera", has_project)
        self._set_toolbar_button_enabled("cut", self._can_cut_selection())
        self._set_toolbar_button_enabled("copy", self._can_copy_selection())
        self._set_toolbar_button_enabled("paste", self._can_paste_selection())
        self._set_toolbar_button_enabled("delete", has_project and has_selection)
        self._set_toolbar_button_enabled("create_complementary", self._can_create_complementary())
        self._set_toolbar_button_enabled("delete_complementary", self._can_delete_complementary())
        self._set_toolbar_button_enabled("create_down", self._can_create_down())
        self._set_toolbar_button_enabled("delete_down", self._can_delete_down())
        self._set_toolbar_button_enabled("create_sticky_end", self._can_create_sticky_end())
        self._set_toolbar_button_enabled("delete_sticky_end", self._can_delete_sticky_end())
        self._set_toolbar_button_enabled("change_backbone_color", self._can_change_backbone_color())
        self._set_toolbar_button_enabled("reset_backbone_color", self._can_reset_backbone_color())
        self._set_toolbar_button_enabled("set_sequence", self._can_set_sequence())
        self._set_toolbar_button_enabled("reset_sequence", self._can_reset_sequence())
        self._set_toolbar_button_enabled("generate_sequence", has_project)
        self._set_toolbar_button_enabled("undo", bool(self.undo_stack))
        self._set_toolbar_button_enabled("redo", bool(self.redo_stack))
        for icon_name in ("move", "move_fixed", "rotate", "rotate_fixed"):
            self._set_toolbar_button_enabled(icon_name, has_project and has_selection)

    def _update_status(self) -> None:
        if not self.project:
            self.status_var.set("Total bases: 0 | Total strands: 0 | Selected: -")
            return
        summary = self.project.summary()
        self.status_var.set(
            f"Total bases: {summary['bases']} | "
            f"Total strands: {summary['strands']} | "
            f"{_selection_status_text(self.project, self.selected_indices, self.primary_selected)}"
        )

    def _redraw_view(self, key: str) -> None:
        state = self.viewports[key]
        canvas = state.canvas
        canvas.delete("all")
        width = max(1, canvas.winfo_width())
        height = max(1, canvas.winfo_height())
        theme = self._theme()

        if not self.project or not self.project.bases:
            self._draw_grid(state, width, height)
            self._draw_corner_badge(state)
            self._draw_legend(state)
            self._draw_axes_indicator(state)
            self._draw_selection_box_overlay(state)
            self._draw_pending_free_strand_overlay(state)
            return

        if state.needs_fit or state.scene_center is None:
            state.scene_center = _project_center(self.project)
        fit_zoom = _fit_camera_zoom(self.project, state.camera, width, height, state.mode, center=state.scene_center)
        if state.needs_fit or state.camera.zoom <= 0.0:
            state.fit_zoom = fit_zoom
            state.camera.zoom = fit_zoom
            state.needs_fit = False
        elif state.fit_zoom <= 0.0:
            state.fit_zoom = fit_zoom

        state.projected = project_points(self.project, state.camera, width, height, mode=state.mode, center=state.scene_center)
        along_edges = self.project.backbone_edges()
        _unused_along_edges, across_edges = self.project.all_edges()
        strand_colors = {strand.strand_id: strand.color for strand in self.project.strands}
        selected_next_along, selected_across = _selected_edge_sets(self.project, self.selected_indices)
        backbone_width = self.display_backbone_width
        basepair_width = self.display_basepair_width
        selected_backbone_width = backbone_width + 2
        radius = _scaled_marker_radius(self.display_sphere_size, state.camera.zoom, state.fit_zoom)
        selected_radius = _scaled_selected_radius(radius)

        self._draw_grid(state, width, height)

        for left, right in across_edges:
            self._draw_edge(
                canvas,
                state.projected,
                left,
                right,
                color="#64748b",
                width=basepair_width,
                dash=_basepair_dash(self.project, left, right),
            )
        for left, right in along_edges:
            self._draw_edge(
                canvas,
                state.projected,
                left,
                right,
                color=_backbone_color_for_base(self.project, left, strand_colors),
                width=backbone_width,
            )
        self._draw_strand_end_markers(state, backbone_width, radius, strand_colors)
        self._draw_sticky_end_overlays(state, radius)
        for left, right in selected_across:
            self._draw_edge(
                canvas,
                state.projected,
                left,
                right,
                color=theme["selection_outline"],
                width=basepair_width,
                dash=_basepair_dash(self.project, left, right),
            )
        for left, right in selected_next_along:
            strand_color = _backbone_color_for_base(self.project, left, strand_colors)
            self._draw_edge(canvas, state.projected, left, right, color=theme["selection_outline"], width=selected_backbone_width)
            self._draw_edge(canvas, state.projected, left, right, color=strand_color, width=backbone_width)

        for index, (x, y, depth) in sorted(state.projected.items(), key=lambda item: item[1][2]):
            base = self.project.bases[index]
            is_selected = index in self.selected_indices
            marker_radius = selected_radius if is_selected else radius
            outline = theme["selection_outline"] if is_selected else theme["base_outline"]
            self._draw_base_marker(
                canvas=canvas,
                base=base,
                x=x,
                y=y,
                radius=marker_radius,
                outline=outline,
                width=2 if is_selected else 1,
            )

        if self.show_strand_labels:
            self._draw_strand_labels(state, radius)
        self._draw_corner_badge(state)
        self._draw_legend(state)
        self._draw_axes_indicator(state)
        self._draw_selection_box_overlay(state)
        self._draw_pending_free_strand_overlay(state)

    def _draw_edge(
        self,
        canvas: tk.Canvas,
        projected: dict[int, tuple[float, float, float]],
        left: int,
        right: int,
        color: str,
        width: int,
        dash: tuple[int, int] | None = None,
    ) -> None:
        x1, y1, _ = projected[left]
        x2, y2, _ = projected[right]
        canvas.create_line(x1, y1, x2, y2, fill=color, width=width, dash=dash)

    def _draw_grid(self, state: ViewportState, width: int, height: int) -> None:
        if not self.show_grid:
            return
        theme = self._theme()
        center = state.scene_center or (0.0, 0.0, 0.0)
        scale = (
            state.camera.zoom
            if state.camera.zoom > 0.0
            else _grid_fit_zoom_for_mode(self.grid_extents, state.camera, width, height, state.mode)
        )
        if scale <= 0.0:
            return
        plane_modes = ("xy", "yz", "xz") if state.mode == "3d" else (state.mode,)
        spacing = self.grid_extents["spacing"]
        grid_color = _blend_hex(theme["grid_line"], theme["canvas_bg"], 0.5)

        for plane_mode in plane_modes:
            axis_x, axis_y = _grid_axes_for_mode(plane_mode)
            for value in _grid_values(self.grid_extents[f"{axis_x}_min"], self.grid_extents[f"{axis_x}_max"], spacing):
                x1, y1 = _project_world_to_canvas(
                    _grid_world_point(plane_mode, value, self.grid_extents[f"{axis_y}_min"]),
                    state.camera,
                    width,
                    height,
                    state.mode,
                    center,
                    scale,
                )
                x2, y2 = _project_world_to_canvas(
                    _grid_world_point(plane_mode, value, self.grid_extents[f"{axis_y}_max"]),
                    state.camera,
                    width,
                    height,
                    state.mode,
                    center,
                    scale,
                )
                color = _blend_hex(AXIS_COLORS[axis_y], theme["canvas_bg"], 0.7) if abs(value) <= spacing * 0.5 else grid_color
                line_width = 1.5 if abs(value) <= spacing * 0.5 else 1.0
                state.canvas.create_line(x1, y1, x2, y2, fill=color, width=line_width)

            for value in _grid_values(self.grid_extents[f"{axis_y}_min"], self.grid_extents[f"{axis_y}_max"], spacing):
                x1, y1 = _project_world_to_canvas(
                    _grid_world_point(plane_mode, self.grid_extents[f"{axis_x}_min"], value),
                    state.camera,
                    width,
                    height,
                    state.mode,
                    center,
                    scale,
                )
                x2, y2 = _project_world_to_canvas(
                    _grid_world_point(plane_mode, self.grid_extents[f"{axis_x}_max"], value),
                    state.camera,
                    width,
                    height,
                    state.mode,
                    center,
                    scale,
                )
                color = _blend_hex(AXIS_COLORS[axis_x], theme["canvas_bg"], 0.7) if abs(value) <= spacing * 0.5 else grid_color
                line_width = 1.5 if abs(value) <= spacing * 0.5 else 1.0
                state.canvas.create_line(x1, y1, x2, y2, fill=color, width=line_width)

    def _draw_base_marker(
        self,
        canvas: tk.Canvas,
        base,
        x: float,
        y: float,
        radius: float,
        outline: str,
        width: int,
    ) -> None:
        fill = BASE_COLORS.get(base.nucleotide, BASE_COLORS[None])
        if str(getattr(base, "molecule", "DNA")).upper() == "RNA":
            self._draw_rna_square(canvas, x, y, radius, fill, outline, width)
            return
        canvas.create_oval(
            x - radius,
            y - radius,
            x + radius,
            y + radius,
            fill=fill,
            outline=outline,
            width=width,
        )

    def _draw_rna_square(
        self,
        canvas: tk.Canvas,
        x: float,
        y: float,
        radius: float,
        fill: str,
        outline: str,
        width: int,
    ) -> None:
        canvas.create_rectangle(
            x - radius,
            y - radius,
            x + radius,
            y + radius,
            fill=fill,
            outline=outline,
            width=width,
        )

    def _draw_strand_labels(self, state: ViewportState, radius: float) -> None:
        if not self.project:
            return
        canvas = state.canvas
        theme = self._theme()
        for strand in self.project.strands:
            if not strand.base_indices:
                continue
            first_index = strand.base_indices[0]
            if first_index not in state.projected:
                continue
            x, y, _depth = state.projected[first_index]
            label_x = x + radius + 2.0
            label_y = y - radius - 2.0
            canvas.create_text(
                label_x,
                label_y,
                anchor="sw",
                fill=theme["label_fg"],
                font=("TkDefaultFont", 10, "bold"),
                text=str(strand.strand_id),
            )

    def _draw_strand_end_markers(
        self,
        state: ViewportState,
        backbone_width: int,
        radius: float,
        strand_colors: dict[int, str],
    ) -> None:
        if not self.project:
            return
        canvas = state.canvas
        extension, arrowshape, _svg_length, _svg_half_width = _strand_end_arrow_style(radius)
        for strand in self.project.strands:
            segment = _strand_end_arrow_segment(self.project, strand, state.projected, extension=extension)
            if segment is None:
                continue
            x1, y1, x2, y2 = segment
            canvas.create_line(
                x1,
                y1,
                x2,
                y2,
                fill=_strand_end_color(self.project, strand, strand_colors),
                width=backbone_width,
                arrow=tk.LAST,
                arrowshape=arrowshape,
            )

    def _draw_sticky_end_overlays(self, state: ViewportState, radius: float) -> None:
        if not self.project or not self.show_sticky_end:
            return
        dash = (12, 8)
        padding = max(4.0, radius + 2.0)
        canvas = state.canvas
        for group in self.project.sticky_end_groups():
            color = str(group["color"])
            for left_index, right_index in _sticky_end_pairs(group):
                if left_index not in state.projected or right_index not in state.projected:
                    continue
                self._draw_edge(
                    canvas,
                    state.projected,
                    left_index,
                    right_index,
                    color=color,
                    width=1,
                    dash=dash,
                )
            for indices in (group["left_indices"], group["right_indices"]):
                bounds = _projected_group_bounds(state.projected, indices, padding)
                if bounds is None:
                    continue
                canvas.create_rectangle(*bounds, outline=color, width=1)

    def _draw_corner_badge(self, state: ViewportState) -> None:
        badge = VIEW_BADGES.get(state.key)
        if not badge:
            return
        height = max(1, state.canvas.winfo_height())
        theme = self._theme()
        state.canvas.create_text(
            18,
            height - 18,
            anchor="sw",
            fill=theme["badge_fg"],
            font=("TkDefaultFont", 14, "bold"),
            text=badge,
        )

    def _draw_legend(self, state: ViewportState) -> None:
        if state.key != "xz":
            return
        x = 18
        y = 16
        for letter, color in LEGEND_ITEMS:
            state.canvas.create_text(
                x,
                y,
                anchor="nw",
                fill=color,
                font=("TkDefaultFont", 11, "bold"),
                text=letter,
            )
            y += 16

    def _draw_axes_indicator(self, state: ViewportState) -> None:
        width = max(1, state.canvas.winfo_width())
        height = max(1, state.canvas.winfo_height())
        origin_x = width - 62
        origin_y = height - 42
        length = 28

        if state.mode == "3d":
            axes = {
                axis_name: rotate_point(*vector, state.camera)
                for axis_name, vector in {
                    "x": (1.0, 0.0, 0.0),
                    "y": (0.0, 1.0, 0.0),
                    "z": (0.0, 0.0, 1.0),
                }.items()
            }
            for axis_name, (ax, ay, _az) in axes.items():
                dx = ax * length
                dy = -ay * length
                self._draw_axis_arrow(state.canvas, origin_x, origin_y, dx, dy, axis_name)
            return

        axis_vectors = {
            "xy": (("x", length, 0), ("y", 0, -length)),
            "yz": (("y", length, 0), ("z", 0, -length)),
            "xz": (("x", length, 0), ("z", 0, -length)),
        }
        for axis_name, dx, dy in axis_vectors[state.mode]:
            self._draw_axis_arrow(state.canvas, origin_x, origin_y, dx, dy, axis_name)

    def _draw_axis_arrow(
        self,
        canvas: tk.Canvas,
        origin_x: float,
        origin_y: float,
        dx: float,
        dy: float,
        axis_name: str,
    ) -> None:
        end_x = origin_x + dx
        end_y = origin_y + dy
        color = AXIS_COLORS[axis_name]
        canvas.create_line(
            origin_x,
            origin_y,
            end_x,
            end_y,
            fill=color,
            width=2,
            arrow=tk.LAST,
            arrowshape=(8, 10, 3),
        )
        label_x = end_x + (6 if dx >= 0 else -6)
        label_y = end_y + (6 if dy >= 0 else -6)
        anchor = "nw"
        if dx < 0 and dy <= 0:
            anchor = "se"
        elif dx < 0:
            anchor = "ne"
        elif dy < 0:
            anchor = "sw"
        canvas.create_text(
            label_x,
            label_y,
            anchor=anchor,
            fill=color,
            font=("TkDefaultFont", 10, "bold"),
            text=axis_name,
        )

    def _draw_selection_box_overlay(self, state: ViewportState) -> None:
        if state.selection_box_start is None or state.selection_box_end is None:
            return
        x0, y0 = state.selection_box_start
        x1, y1 = state.selection_box_end
        theme = self._theme()
        if self.selection_mode == "create_helix":
            state.canvas.create_line(
                x0,
                y0,
                x1,
                y1,
                fill=theme["selection_box"],
                width=2,
                dash=(6, 3),
            )
            return
        state.canvas.create_rectangle(
            x0,
            y0,
            x1,
            y1,
            outline=theme["selection_box"],
            width=1,
            dash=(4, 2),
        )

    def _draw_pending_free_strand_overlay(self, state: ViewportState) -> None:
        if self.selection_mode != "create_base" or self.pending_free_strand_start is None:
            return
        width = max(1, state.canvas.winfo_width())
        height = max(1, state.canvas.winfo_height())
        center = state.scene_center or (0.0, 0.0, 0.0)
        scale = _view_scale(self.project, self.grid_extents, state, width, height, center)
        start_x, start_y = _project_world_to_canvas(
            self.pending_free_strand_start.point,
            state.camera,
            width,
            height,
            state.mode,
            center,
            scale,
        )
        theme = self._theme()
        if self.pending_free_strand_start.base_index is None:
            state.canvas.create_oval(
                start_x - 4,
                start_y - 4,
                start_x + 4,
                start_y + 4,
                outline=theme["selection_box"],
                width=2,
            )
        if self.pending_free_strand_hover_view == state.key and self.pending_free_strand_hover_canvas is not None:
            end_x, end_y = self.pending_free_strand_hover_canvas
            state.canvas.create_line(
                start_x,
                start_y,
                end_x,
                end_y,
                fill=theme["selection_box"],
                width=2,
                dash=(6, 3),
            )

    def _begin_drag(self, key: str, button: int, event: tk.Event) -> None:
        state = self.viewports[key]
        state.canvas.focus_set()
        state.drag_anchor = (event.x, event.y)
        state.drag_camera = _copy_camera(state.camera)
        state.drag_button = button
        state.drag_moved = False
        state.drag_positions = None
        state.drag_center = None
        state.drag_projected_center = None
        state.drag_snapshot = None
        if button == 1 and self.selection_mode in {"box", "create_helix"}:
            state.selection_box_start = (event.x, event.y)
            state.selection_box_end = (event.x, event.y)
        else:
            state.selection_box_start = None
            state.selection_box_end = None
        if button == 1 and self.selection_mode in TRANSFORM_MODE_NAMES and self.project and self.selected_indices:
            state.drag_snapshot = self._capture_history_snapshot()
            state.drag_positions = {
                index: self.project.bases[index].position
                for index in self.selected_indices
            }
            state.drag_center = _selection_center(self.project, self.selected_indices)
            state.drag_projected_center = _projected_selection_center(state.projected, self.selected_indices)

    def _drag_view(self, key: str, button: int, event: tk.Event) -> None:
        state = self.viewports[key]
        if state.drag_anchor is None or state.drag_camera is None or state.drag_button != button:
            return
        dx = event.x - state.drag_anchor[0]
        dy = event.y - state.drag_anchor[1]
        if button == 1 and self.selection_mode in {"box", "create_helix"}:
            if abs(dx) + abs(dy) >= 4:
                state.drag_moved = True
            if self.selection_mode == "create_helix":
                state.selection_box_end = self._create_helix_drag_endpoint(key, state, event)
            else:
                state.selection_box_end = (event.x, event.y)
            self._redraw_view(key)
            return
        if button == 1 and self.selection_mode in TRANSFORM_MODE_NAMES and self.project and state.drag_positions:
            if abs(dx) + abs(dy) >= 2:
                state.drag_moved = True
            constrain_axis = _has_shift_modifier(event.state)
            if self.selection_mode == "move":
                self._drag_move_selection(state, dx, dy, constrain_axis)
            else:
                self._drag_rotate_selection(state, dx, dy, constrain_axis)
            return
        if button != 1 and abs(dx) + abs(dy) >= 4:
            state.drag_moved = True

        if button == ROTATE_DRAG_BUTTON and state.mode == "3d":
            state.camera.orientation = _apply_trackball_rotation(state.drag_camera, dx, dy)
            self._redraw_view(key)
            return

        if button == PAN_DRAG_BUTTON:
            state.camera.pan_x = state.drag_camera.pan_x + dx
            state.camera.pan_y = state.drag_camera.pan_y + dy
            self._redraw_view(key)

    def _end_drag(self, key: str, button: int, event: tk.Event) -> None:
        state = self.viewports[key]
        if state.drag_button == button and button == 1:
            if self.selection_mode == "box" and state.drag_moved:
                targets = self._selection_box_targets(key)
                state.selection_box_start = None
                state.selection_box_end = None
                self._apply_selection_targets(targets, primary_index=min(targets) if targets else None, shift=_has_shift_modifier(event.state))
            elif self.selection_mode == "create_helix":
                start_canvas = state.selection_box_start or (event.x, event.y)
                end_canvas = state.selection_box_end or (event.x, event.y)
                state.selection_box_start = None
                state.selection_box_end = None
                self._create_helix_from_interaction(key, start_canvas, end_canvas, state.drag_moved)
            elif self.selection_mode in TRANSFORM_MODE_NAMES:
                state.selection_box_start = None
                state.selection_box_end = None
                if state.drag_moved:
                    self._commit_history_snapshot(state.drag_snapshot)
                    self._refresh_all()
                else:
                    self._redraw_view(key)
            elif self.selection_mode == "create_base" and not state.drag_moved:
                state.selection_box_start = None
                state.selection_box_end = None
                self._handle_create_free_strand_click(key, event)
            elif not state.drag_moved:
                state.selection_box_start = None
                state.selection_box_end = None
                self._select_nearest(key, event)
            else:
                state.selection_box_start = None
                state.selection_box_end = None
                self._redraw_view(key)
        state.drag_anchor = None
        state.drag_camera = None
        state.drag_button = None
        state.drag_moved = False
        state.drag_positions = None
        state.drag_center = None
        state.drag_projected_center = None
        state.drag_snapshot = None

    def _drag_move_selection(self, state: ViewportState, dx: float, dy: float, constrain_axis: bool = False) -> None:
        if not self.project or state.drag_camera is None or state.drag_positions is None:
            return
        scale = state.drag_camera.zoom if state.drag_camera.zoom > 0.0 else max(state.fit_zoom, 1.0)
        translation = _drag_translation_vector(state.mode, state.drag_camera, dx, dy, scale, constrain_axis=constrain_axis)
        _apply_translation_to_project(self.project, state.drag_positions, translation)
        self._refresh_views_only()

    def _drag_rotate_selection(self, state: ViewportState, dx: float, dy: float, constrain_axis: bool = False) -> None:
        if not self.project or state.drag_positions is None or state.drag_center is None or state.drag_camera is None:
            return
        if state.mode == "3d":
            rotation = _drag_rotation_matrix_3d(state.drag_camera, dx, dy, constrain_axis=constrain_axis)
        else:
            rotation = _drag_rotation_matrix_planar(state.mode, dx, dy, constrain_axis=constrain_axis)
        _apply_rotation_to_project(self.project, state.drag_positions, state.drag_center, rotation)
        self._refresh_views_only()

    def _select_nearest(self, key: str, event: tk.Event) -> None:
        if not self.project:
            return
        state = self.viewports[key]
        best_index = _nearest_projected_index(state.projected, event.x, event.y)
        if best_index is None:
            if _clear_selection_on_empty_click(self.selection_mode):
                self.selected_indices.clear()
                self.primary_selected = None
                self._refresh_all()
            return
        targets = _selection_targets_for_mode(self.project, best_index, self.selection_mode)
        self._apply_selection_targets(targets, primary_index=best_index, shift=_has_shift_modifier(event.state))

    def _queue_base_editor_from_click(self, key: str, event: tk.Event) -> str:
        if self.selection_mode in {"create_helix", "create_base"} or not self.project:
            return "break"
        state = self.viewports[key]
        best_index = _nearest_projected_index(state.projected, event.x, event.y)
        if best_index is None:
            return "break"
        self.root.after_idle(lambda selected_index=best_index: self._open_base_editor(selected_index))
        return "break"

    def _open_base_editor(self, index: int) -> None:
        if not self.project or index not in self.project.bases:
            return
        base = self.project.bases[index]
        self.selected_indices = {index}
        self.primary_selected = index
        self._refresh_all()

        options = _base_edit_options(base.molecule)
        current_value = base.nucleotide if base.nucleotide in {value for _label, value in options} else "N"

        dialog = tk.Toplevel(self.root)
        dialog.withdraw()
        dialog.title("Change Base")
        dialog.configure(bg="#e5e5e5")
        dialog.resizable(False, False)
        dialog.transient(self.root)

        value_var = tk.StringVar(value=current_value)
        changed = {"value": False}

        body = tk.Frame(dialog, bg="#e5e5e5", padx=12, pady=12)
        body.pack(fill=tk.BOTH, expand=True)

        options_frame = tk.Frame(body, bg="#e5e5e5")
        options_frame.grid(row=0, column=0, sticky="nw", padx=(0, 18))

        initial_focus = None
        for row, (label, value) in enumerate(options):
            radio = tk.Radiobutton(
                options_frame,
                text=label,
                value=value,
                variable=value_var,
                bg="#e5e5e5",
                activebackground="#e5e5e5",
                highlightthickness=0,
                anchor="w",
                justify=tk.LEFT,
                selectcolor="#ffffff",
            )
            radio.grid(row=row, column=0, sticky="w", pady=2)
            if value == current_value and initial_focus is None:
                initial_focus = radio

        buttons = tk.Frame(body, bg="#e5e5e5")
        buttons.grid(row=0, column=1, sticky="ne")

        def confirm() -> None:
            new_value = value_var.get()
            new_nucleotide = None if new_value == "N" else new_value
            if base.nucleotide != new_nucleotide:
                snapshot = self._capture_history_snapshot()
                self.project.set_nucleotide(index, new_nucleotide)
                self._commit_history_snapshot(snapshot)
                changed["value"] = True
            dialog.destroy()

        ok_button = tk.Button(buttons, text="OK", width=9, command=confirm, highlightthickness=0)
        ok_button.grid(row=0, column=0, sticky="ew", pady=(0, 8))
        cancel_button = tk.Button(buttons, text="Cancel", width=9, command=dialog.destroy, highlightthickness=0)
        cancel_button.grid(row=1, column=0, sticky="ew")

        dialog.bind("<Return>", lambda _event: confirm())
        dialog.bind("<Escape>", lambda _event: dialog.destroy())
        dialog.protocol("WM_DELETE_WINDOW", dialog.destroy)

        self.root.update_idletasks()
        dialog.update_idletasks()
        x = self.root.winfo_rootx() + max((self.root.winfo_width() - dialog.winfo_width()) // 2, 0)
        y = self.root.winfo_rooty() + max((self.root.winfo_height() - dialog.winfo_height()) // 2, 0)
        dialog.geometry(f"+{x}+{y}")
        dialog.deiconify()
        dialog.lift()
        dialog.grab_set()
        if initial_focus is not None:
            initial_focus.focus_set()
        else:
            dialog.focus_set()
        dialog.wait_window()
        if changed["value"]:
            self._refresh_all()

    def _selection_box_targets(self, key: str) -> set[int]:
        state = self.viewports[key]
        if state.selection_box_start is None or state.selection_box_end is None:
            return set()
        return _projected_indices_in_box(state.projected, state.selection_box_start, state.selection_box_end)

    def _apply_selection_targets(self, targets: set[int], primary_index: int | None, shift: bool) -> None:
        if not self.project:
            return
        valid_targets = {index for index in targets if index in self.project.bases}
        if shift:
            if not valid_targets:
                return
            if valid_targets.issubset(self.selected_indices):
                self.selected_indices.difference_update(valid_targets)
                if self.primary_selected in valid_targets:
                    self.primary_selected = min(self.selected_indices) if self.selected_indices else None
            else:
                self.selected_indices.update(valid_targets)
                if primary_index in valid_targets:
                    self.primary_selected = primary_index
                else:
                    self.primary_selected = min(valid_targets)
        else:
            self.selected_indices = set(valid_targets)
            if primary_index in self.selected_indices:
                self.primary_selected = primary_index
            else:
                self.primary_selected = min(self.selected_indices) if self.selected_indices else None
        self._refresh_all()

    def _on_wheel(self, key: str, event: tk.Event) -> None:
        if self._zoom_blocked(key):
            return "break"
        delta = 0.1 if event.delta > 0 else -0.1
        self._zoom_view(key, delta)
        return "break"

    def _on_scroll_zoom(self, key: str, delta: float) -> str | None:
        if self._zoom_blocked(key):
            return "break"
        self._zoom_view(key, delta)
        return "break"

    def _zoom_view(self, key: str, delta: float) -> None:
        state = self.viewports[key]
        width = max(1, state.canvas.winfo_width())
        height = max(1, state.canvas.winfo_height())
        if state.needs_fit or state.camera.zoom <= 0.0:
            if self.project and self.project.bases:
                if state.scene_center is None:
                    state.scene_center = _project_center(self.project)
                state.fit_zoom = _fit_camera_zoom(self.project, state.camera, width, height, state.mode, center=state.scene_center)
            else:
                state.scene_center = (0.0, 0.0, 0.0)
                state.fit_zoom = _grid_fit_zoom_for_mode(self.grid_extents, state.camera, width, height, state.mode)
            state.camera.zoom = state.fit_zoom
            state.needs_fit = False
        zoom_factor = 1.2 ** (delta / 0.1)
        state.camera.zoom = max(ABSOLUTE_ZOOM_MIN, min(ABSOLUTE_ZOOM_MAX, state.camera.zoom * zoom_factor))
        self._redraw_view(key)
        self._update_status()

    def _zoom_blocked(self, key: str) -> bool:
        state = self.viewports[key]
        return state.drag_button is not None

    def _on_keypress(self, event: tk.Event) -> str | None:
        mode = MODE_KEYS.get(event.keysym)
        if mode is not None and not _has_shortcut_modifier(event.state):
            self._set_selection_mode(mode)
            return "break"
        if not self.project:
            return None
        if _has_shortcut_modifier(event.state):
            return None
        if not self.selected_indices:
            return None

        key = event.keysym.upper()
        if key in {"A", "C", "G", "T", "U"}:
            snapshot = self._capture_history_snapshot()
            self._set_selection_nucleotide(key)
            self._commit_history_snapshot(snapshot)
            self._refresh_all()
            return "break"
        if key in {"BACKSPACE", "DELETE"}:
            self._delete_selection()
            return "break"
        return None

    def _set_selection_nucleotide(self, nucleotide: str) -> None:
        if not self.project:
            return
        selected = set(self.selected_indices)
        for index in selected:
            self.project.set_nucleotide(index, nucleotide)

    def _select_all(self) -> None:
        if not self.project:
            return
        self.selected_indices = set(self.project.bases)
        self.primary_selected = min(self.selected_indices) if self.selected_indices else None
        self._refresh_all()

    def _select_all_event(self, _event: tk.Event) -> str:
        self._select_all()
        return "break"

    def _new_project(self) -> None:
        self.project = None
        self.constraints_values = _default_constraints_values()
        self.custom_constraints_values = _blank_custom_constraints_values()
        self.selected_indices.clear()
        self.primary_selected = None
        self._clear_pending_free_strand()
        self._clear_history()
        self._reset_viewports()
        self._refresh_all()

    def _new_project_event(self, _event: tk.Event) -> str:
        self._new_project()
        return "break"

    def _open_file_event(self, _event: tk.Event) -> str:
        self.open_file()
        return "break"

    def _save_project_event(self, _event: tk.Event) -> str:
        self.save_project()
        return "break"

    def _save_project_as_event(self, _event: tk.Event) -> str:
        self.save_json()
        return "break"

    def _cut_event(self, _event: tk.Event) -> str:
        self._cut_selection()
        return "break"

    def _copy_event(self, _event: tk.Event) -> str:
        self._copy_selection()
        return "break"

    def _paste_event(self, _event: tk.Event) -> str:
        self._paste_selection()
        return "break"

    def _undo_event(self, _event: tk.Event) -> str:
        self.undo()
        return "break"

    def _redo_event(self, _event: tk.Event) -> str:
        self.redo()
        return "break"

    def _delete_selection(self) -> None:
        if not self.project or not self.selected_indices:
            return
        snapshot = self._capture_history_snapshot()
        self.project.delete_bases(self.selected_indices)
        self.selected_indices.clear()
        self.primary_selected = None
        self._commit_history_snapshot(snapshot)
        self._refresh_all()

    def open_file(self) -> None:
        path = filedialog.askopenfilename(
            title="Open Tiamat data",
            filetypes=[
                ("Tiamat binary archive", "*.dna"),
                ("Tiamat JSON", "*.dnajson *.json"),
                ("Nucleic acid PDB", "*.pdb"),
                ("All files", "*.*"),
            ],
        )
        if not path:
            return
        try:
            self._open_project_path(path)
        except Exception as exc:
            messagebox.showerror("Open failed", str(exc))

    def save_project(self) -> None:
        if not self.project:
            return
        current_path = self._project_save_path()
        if current_path is None:
            self.save_json()
            return
        self._write_project_json(current_path)

    def save_json(self) -> None:
        if not self.project:
            return
        path = filedialog.asksaveasfilename(
            title="Save As",
            defaultextension=".dnajson",
            filetypes=[("Tiamat JSON", "*.dnajson"), ("JSON", "*.json")],
        )
        if not path:
            return
        self._write_project_json(path)

    def export_sequences(self) -> None:
        if not self.project:
            return
        path = filedialog.asksaveasfilename(
            title="Export sequences",
            defaultextension=".txt",
            filetypes=[("Text", "*.txt")],
        )
        if not path:
            return
        Path(path).write_text(self.project.export_sequences_text(), encoding="utf-8")

    def render_project(self) -> None:
        if not self.project:
            return
        options = self._prompt_render_dialog()
        if options is None:
            return
        destination = options["destination"]
        if destination == "Clipboard":
            try:
                self._copy_render_png_to_clipboard(options["view"])
            except Exception as exc:
                messagebox.showerror("Render", str(exc), parent=self.root)
                return
            messagebox.showinfo("Render", f"{options['view']} PNG copied to the clipboard.", parent=self.root)
            return
        svg, default_name = self._render_svg_for_view(options["view"])
        path = filedialog.asksaveasfilename(
            title="Render",
            initialfile=default_name,
            defaultextension=".svg",
            filetypes=[("SVG", "*.svg")],
        )
        if not path:
            return
        Path(path).write_text(svg, encoding="utf-8")

    def _copy_render_png_to_clipboard(self, view_label: str) -> None:
        if sys.platform != "darwin":
            raise RuntimeError("Clipboard image rendering is currently supported on macOS only.")
        width, height = self._render_base_dimensions(view_label)
        scale = _clipboard_render_scale(width, height, max_dimension=CLIPBOARD_RENDER_MAX_DIMENSION)
        svg, default_name = self._render_svg_for_view(view_label, scale=scale)
        with tempfile.TemporaryDirectory(prefix="tiamat_py_render_") as tmpdir:
            temp_dir = Path(tmpdir)
            svg_path = temp_dir / default_name
            png_path = svg_path.with_suffix(".png")
            svg_path.write_text(svg, encoding="utf-8")
            try:
                subprocess.run(
                    ["sips", "-s", "format", "png", str(svg_path), "--out", str(png_path)],
                    check=True,
                    capture_output=True,
                    text=True,
                )
            except OSError as exc:
                raise RuntimeError(f"Failed to start PNG conversion for the clipboard. {exc}") from exc
            except subprocess.CalledProcessError as exc:
                details = exc.stderr.strip() or exc.stdout.strip() or "Unknown sips error."
                raise RuntimeError(f"Failed to render PNG for the clipboard. {details}") from exc
            swift_script = """
import AppKit

let path = CommandLine.arguments[1]
guard let image = NSImage(contentsOfFile: path) else {
    fputs("Could not load rendered PNG.\\n", stderr)
    exit(1)
}
let pasteboard = NSPasteboard.general
pasteboard.clearContents()
if !pasteboard.writeObjects([image]) {
    fputs("Could not write PNG to the clipboard.\\n", stderr)
    exit(2)
}
"""
            try:
                subprocess.run(
                    ["swift", "-", str(png_path)],
                    input=swift_script,
                    check=True,
                    capture_output=True,
                    text=True,
                )
            except OSError as exc:
                raise RuntimeError(f"Failed to start clipboard image transfer. {exc}") from exc
            except subprocess.CalledProcessError as exc:
                details = exc.stderr.strip() or exc.stdout.strip() or "Unknown Swift clipboard error."
                raise RuntimeError(f"Failed to copy the rendered PNG to the clipboard. {details}") from exc

    def _render_base_dimensions(self, view_label: str) -> tuple[int, int]:
        if not self.project:
            raise ValueError("No project loaded.")
        viewport_key = dict(RENDER_VIEW_OPTIONS).get(view_label, "view3d")
        state = self.viewports[viewport_key]
        width = state.canvas.winfo_width()
        height = state.canvas.winfo_height()
        if width <= 1 or height <= 1:
            return (1200, 900)
        return (width, height)

    def _render_svg_for_view(self, view_label: str, scale: float = 1.0) -> tuple[str, str]:
        if not self.project:
            raise ValueError("No project loaded.")
        viewport_key = dict(RENDER_VIEW_OPTIONS).get(view_label, "view3d")
        state = self.viewports[viewport_key]
        base_width, base_height = self._render_base_dimensions(view_label)
        scale = max(scale, 1.0)
        width = max(1, round(base_width * scale))
        height = max(1, round(base_height * scale))
        camera = _copy_camera(state.camera)
        camera.zoom *= scale
        camera.pan_x *= scale
        camera.pan_y *= scale
        fit_zoom = state.fit_zoom
        if fit_zoom <= 0.0:
            fit_zoom = _fit_camera_zoom(self.project, state.camera, base_width, base_height, state.mode, center=state.scene_center)
        svg = render_svg(
            self.project,
            camera,
            width,
            height,
            self.primary_selected,
            selected_indices=self.selected_indices,
            sphere_size=self.display_sphere_size * scale,
            backbone_width=max(1, round(self.display_backbone_width * scale)),
            basepair_width=max(1, round(self.display_basepair_width * scale)),
            fit_zoom=fit_zoom * scale,
            center=state.scene_center,
            color_scheme=self.color_scheme,
            show_sticky_end=self.show_sticky_end,
            mode=state.mode,
        )
        source_path = self.project.metadata.get("source_path")
        stem = Path(str(source_path)).stem if source_path else "render"
        default_name = f"{stem}_{view_label.lower()}.svg"
        return svg, default_name

    def _prompt_render_dialog(self) -> dict[str, str] | None:
        dialog = tk.Toplevel(self.root)
        dialog.withdraw()
        dialog.title("Render")
        dialog.configure(bg="#e5e5e5")
        dialog.resizable(False, False)
        dialog.transient(self.root)

        view_var = tk.StringVar(value="Perspective")
        destination_var = tk.StringVar(value="File")
        result: dict[str, str] | None = None

        body = tk.Frame(dialog, bg="#e5e5e5", padx=12, pady=12)
        body.pack(fill=tk.BOTH, expand=True)

        view_frame = tk.Frame(body, bg="#e5e5e5")
        view_frame.grid(row=0, column=0, sticky="nw", padx=(0, 18))
        tk.Label(view_frame, text="View", bg="#e5e5e5", anchor="w").grid(row=0, column=0, sticky="w")
        initial_focus = None
        for row, (label, _viewport_key) in enumerate(RENDER_VIEW_OPTIONS, start=1):
            radio = tk.Radiobutton(
                view_frame,
                text=label,
                value=label,
                variable=view_var,
                bg="#e5e5e5",
                activebackground="#e5e5e5",
                highlightthickness=0,
                anchor="w",
                justify=tk.LEFT,
                selectcolor="#ffffff",
            )
            radio.grid(row=row, column=0, sticky="w", pady=2)
            if initial_focus is None:
                initial_focus = radio

        destination_frame = tk.Frame(body, bg="#e5e5e5")
        destination_frame.grid(row=0, column=1, sticky="nw", padx=(0, 18))
        tk.Label(destination_frame, text="Destination", bg="#e5e5e5", anchor="w").grid(row=0, column=0, sticky="w")
        for row, label in enumerate(RENDER_DESTINATION_OPTIONS, start=1):
            tk.Radiobutton(
                destination_frame,
                text=label,
                value=label,
                variable=destination_var,
                bg="#e5e5e5",
                activebackground="#e5e5e5",
                highlightthickness=0,
                anchor="w",
                justify=tk.LEFT,
                selectcolor="#ffffff",
            ).grid(row=row, column=0, sticky="w", pady=2)

        buttons = tk.Frame(body, bg="#e5e5e5")
        buttons.grid(row=0, column=2, sticky="ne")

        def confirm() -> None:
            nonlocal result
            result = {
                "view": view_var.get(),
                "destination": destination_var.get(),
            }
            dialog.destroy()

        ok_button = tk.Button(buttons, text="OK", width=9, command=confirm, highlightthickness=0)
        ok_button.grid(row=0, column=0, sticky="ew", pady=(0, 8))
        cancel_button = tk.Button(buttons, text="Cancel", width=9, command=dialog.destroy, highlightthickness=0)
        cancel_button.grid(row=1, column=0, sticky="ew")

        dialog.bind("<Return>", lambda _event: confirm())
        dialog.bind("<Escape>", lambda _event: dialog.destroy())
        dialog.protocol("WM_DELETE_WINDOW", dialog.destroy)

        self.root.update_idletasks()
        dialog.update_idletasks()
        x = self.root.winfo_rootx() + max((self.root.winfo_width() - dialog.winfo_width()) // 2, 0)
        y = self.root.winfo_rooty() + max((self.root.winfo_height() - dialog.winfo_height()) // 2, 0)
        dialog.geometry(f"+{x}+{y}")
        dialog.deiconify()
        dialog.lift()
        dialog.grab_set()
        if initial_focus is not None:
            initial_focus.focus_set()
        dialog.wait_window()
        return result

    def generate_sequences(self) -> None:
        if not self.project:
            return
        options = self._prompt_generate_sequences_dialog()
        if options is None:
            return
        snapshot = self._capture_history_snapshot()
        try:
            self.project.generate_sequences(
                unique_sequence_limit=options["unique_sequence_limit"],
                repetition_limit=options["repetition_limit"],
                gc_percentage=options["gc_percentage"],
                timeout=options["timeout"],
            )
        except Exception as exc:
            messagebox.showerror("Generate Sequence", str(exc), parent=self.root)
            return
        self._commit_history_snapshot(snapshot)
        self._refresh_all()

    def reset_bases(self) -> None:
        if not self.project:
            return
        snapshot = self._capture_history_snapshot()
        targets = sorted(self.selected_indices) if self.selected_indices else None
        self.project.reset_bases(targets)
        self._commit_history_snapshot(snapshot)
        self._refresh_all()

    def _reset_viewports(self) -> None:
        for state in self.viewports.values():
            state.camera = _default_camera(state.mode)
            state.projected.clear()
            state.drag_anchor = None
            state.drag_camera = None
            state.drag_button = None
            state.drag_moved = False
            state.drag_positions = None
            state.drag_center = None
            state.drag_projected_center = None
            state.drag_snapshot = None
            state.selection_box_start = None
            state.selection_box_end = None
            _mark_viewport_for_fit(state)

    def _prompt_generate_sequences_dialog(self) -> dict[str, float | int] | None:
        dialog = tk.Toplevel(self.root)
        dialog.withdraw()
        dialog.title("Generate Sequences")
        dialog.configure(bg="#e5e5e5")
        dialog.resizable(False, False)
        dialog.transient(self.root)

        values = {
            "unique_sequence_limit": tk.StringVar(value="10"),
            "repetition_limit": tk.StringVar(value="6"),
            "gc_percentage": tk.StringVar(value="0.5"),
            "timeout": tk.StringVar(value="16"),
        }
        error_var = tk.StringVar(value="")
        result: dict[str, float | int] | None = None

        body = tk.Frame(dialog, bg="#e5e5e5", padx=12, pady=12)
        body.pack(fill=tk.BOTH, expand=True)

        fields = (
            ("Unique Sequence Limit:", "unique_sequence_limit"),
            ("Repetition Limit:", "repetition_limit"),
            ("G/C Percentage:", "gc_percentage"),
            ("Operation Timeout:", "timeout"),
        )
        first_entry = None
        for row, (label_text, key) in enumerate(fields):
            tk.Label(body, text=label_text, bg="#e5e5e5", anchor="w").grid(row=row, column=0, sticky="w", pady=3)
            entry = tk.Entry(body, textvariable=values[key], width=12)
            entry.grid(row=row, column=1, sticky="ew", padx=(10, 0), pady=3)
            if first_entry is None:
                first_entry = entry

        buttons = tk.Frame(body, bg="#e5e5e5")
        buttons.grid(row=0, column=2, rowspan=len(fields), sticky="ne", padx=(12, 0))
        ok_button = tk.Button(buttons, text="OK", width=9, highlightthickness=0)
        ok_button.grid(row=0, column=0, sticky="ew", pady=(0, 8))
        cancel_button = tk.Button(buttons, text="Cancel", width=9, command=dialog.destroy, highlightthickness=0)
        cancel_button.grid(row=1, column=0, sticky="ew")

        tk.Label(
            body,
            text="G/C Percentage uses a fixed +/- 0.05 tolerance.",
            bg="#e5e5e5",
            anchor="w",
            justify=tk.LEFT,
        ).grid(row=len(fields), column=0, columnspan=2, sticky="w", pady=(8, 0))
        error_label = tk.Label(
            body,
            textvariable=error_var,
            bg="#e5e5e5",
            fg="#b91c1c",
            anchor="w",
            justify=tk.LEFT,
        )
        error_label.grid(row=len(fields) + 1, column=0, columnspan=3, sticky="w", pady=(4, 0))

        def parse() -> dict[str, float | int]:
            unique_sequence_limit = int(values["unique_sequence_limit"].get().strip())
            repetition_limit = int(values["repetition_limit"].get().strip())
            gc_percentage = float(values["gc_percentage"].get().strip())
            timeout = float(values["timeout"].get().strip())
            if unique_sequence_limit < 4:
                raise ValueError("Unique Sequence Limit must be at least 4.")
            if repetition_limit < 2:
                raise ValueError("Repetition Limit must be at least 2.")
            if not 0.0 <= gc_percentage <= 1.0:
                raise ValueError("G/C Percentage must be between 0 and 1.")
            if timeout <= 0.0:
                raise ValueError("Operation Timeout must be greater than 0.")
            return {
                "unique_sequence_limit": unique_sequence_limit,
                "repetition_limit": repetition_limit,
                "gc_percentage": gc_percentage,
                "timeout": timeout,
            }

        def confirm() -> None:
            nonlocal result
            try:
                result = parse()
            except ValueError as exc:
                error_var.set(str(exc))
                return
            dialog.destroy()

        ok_button.configure(command=confirm)
        dialog.bind("<Return>", lambda _event: confirm())
        dialog.bind("<Escape>", lambda _event: dialog.destroy())
        dialog.protocol("WM_DELETE_WINDOW", dialog.destroy)

        self.root.update_idletasks()
        dialog.update_idletasks()
        x = self.root.winfo_rootx() + max((self.root.winfo_width() - dialog.winfo_width()) // 2, 0)
        y = self.root.winfo_rooty() + max((self.root.winfo_height() - dialog.winfo_height()) // 2, 0)
        dialog.geometry(f"+{x}+{y}")
        dialog.deiconify()
        dialog.lift()
        dialog.grab_set()
        if first_entry is not None:
            first_entry.focus_set()
            first_entry.selection_range(0, tk.END)
        dialog.wait_window()
        return result


def load_project(path: str | Path) -> TiamatProject:
    file_path = Path(path)
    suffix = file_path.suffix.lower()
    if suffix in {".json", ".dnajson"}:
        return load_json_project(file_path)
    if suffix == ".dna":
        return load_dna_project(file_path)
    if suffix == ".pdb":
        return load_pdb_sugar_center_project(file_path)
    raise ValueError(f"unsupported input type: {file_path.suffix}; use .dna, .dnajson, or .pdb")


def _snapshot_history_state(
    project: TiamatProject | None,
    selected_indices: set[int],
    primary_selected: int | None,
) -> dict[str, object]:
    project_payload = json.loads(project.to_json()) if project is not None else None
    return {
        "project": project_payload,
        "selected_indices": sorted(selected_indices),
        "primary_selected": primary_selected,
    }


def _restore_history_state(
    snapshot: dict[str, object],
) -> tuple[TiamatProject | None, set[int], int | None]:
    payload = snapshot.get("project")
    project = TiamatProject.from_dict(payload) if isinstance(payload, dict) else None
    selected_payload = snapshot.get("selected_indices")
    selected_indices = {
        int(value)
        for value in selected_payload
    } if isinstance(selected_payload, list) else set()
    primary_selected = snapshot.get("primary_selected")
    primary_index = int(primary_selected) if primary_selected is not None else None
    return project, selected_indices, primary_index


def _selection_clipboard_payload(project: TiamatProject, selected_indices: set[int]) -> dict[str, object] | None:
    indices = sorted(index for index in selected_indices if index in project.bases)
    if not indices:
        return None
    selected = set(indices)
    bases_payload: list[dict[str, object]] = []
    for index in indices:
        base = project.bases[index]
        payload = base.to_dict()
        payload["up"] = base.up if base.up in selected else None
        payload["down"] = base.down if base.down in selected else None
        payload["across"] = base.across if base.across in selected else None
        if payload["down"] is None:
            payload["backbone_color"] = None
        bases_payload.append(payload)
    return {"bases": bases_payload}


def _clipboard_paste_offset(paste_count: int) -> tuple[float, float, float]:
    _unused = paste_count
    return (0.0, 0.0, 0.0)


def _paste_structure_clipboard(
    project: TiamatProject | None,
    clipboard_payload: dict[str, object] | None,
    offset: tuple[float, float, float],
) -> tuple[TiamatProject, list[int]]:
    if not clipboard_payload:
        return project if project is not None else TiamatProject(bases={}), []
    payload_bases = clipboard_payload.get("bases")
    if not isinstance(payload_bases, list) or not payload_bases:
        return project if project is not None else TiamatProject(bases={}), []

    target_project = project if project is not None else TiamatProject(bases={})
    next_index = max(target_project.bases, default=0) + 1
    remap: dict[int, int] = {}
    for item in payload_bases:
        if isinstance(item, dict):
            remap[int(item["index"])] = next_index
            next_index += 1

    new_indices: list[int] = []
    dx, dy, dz = offset
    for item in payload_bases:
        if not isinstance(item, dict):
            continue
        source_index = int(item["index"])
        mapped_index = remap[source_index]
        payload = dict(item)
        payload["index"] = mapped_index
        payload["object_id"] = str(mapped_index)
        payload["x"] = float(item["x"]) + dx
        payload["y"] = float(item["y"]) + dy
        payload["z"] = float(item["z"]) + dz
        payload["up"] = remap.get(int(item["up"])) if item.get("up") is not None else None
        payload["down"] = remap.get(int(item["down"])) if item.get("down") is not None else None
        payload["across"] = remap.get(int(item["across"])) if item.get("across") is not None else None
        if payload["down"] is None:
            payload["backbone_color"] = None
        target_project.bases[mapped_index] = Base.from_dict(payload)
        new_indices.append(mapped_index)

    target_project.refresh_strands()
    return target_project, new_indices


def project_points(
    project: TiamatProject,
    camera: Camera,
    width: int,
    height: int,
    mode: str = "3d",
    center: tuple[float, float, float] | None = None,
) -> dict[int, tuple[float, float, float]]:
    centered = _centered_coordinates(project, center=center)
    scale = camera.zoom if camera.zoom > 0.0 else _fit_camera_zoom(project, camera, width, height, mode, center=center)

    projected: dict[int, tuple[float, float, float]] = {}
    for index, x, y, z in centered:
        rx, ry, rz = _project_vector(x, y, z, camera, mode)
        projected[index] = (
            width / 2.0 + camera.pan_x + rx * scale,
            height / 2.0 + camera.pan_y - ry * scale,
            rz,
        )
    return projected


def rotate_point(x: float, y: float, z: float, camera: Camera) -> tuple[float, float, float]:
    if camera.orientation is not None:
        return _apply_orientation(camera.orientation, x, y, z)
    cx, sx = cos(radians(camera.rot_x)), sin(radians(camera.rot_x))
    cy, sy = cos(radians(camera.rot_y)), sin(radians(camera.rot_y))
    cz, sz = cos(radians(camera.rot_z)), sin(radians(camera.rot_z))

    y, z = y * cx - z * sx, y * sx + z * cx
    x, z = x * cy + z * sy, -x * sy + z * cy
    x, y = x * cz - y * sz, x * sz + y * cz
    return x, y, z


def render_svg(
    project: TiamatProject,
    camera: Camera,
    width: int = 1200,
    height: int = 900,
    selected_index: int | None = None,
    selected_indices: set[int] | None = None,
    sphere_size: float = 4.0,
    backbone_width: int = 2,
    basepair_width: int = 1,
    fit_zoom: float | None = None,
    center: tuple[float, float, float] | None = None,
    color_scheme: str = "dark",
    show_sticky_end: bool = True,
    mode: str = "3d",
) -> str:
    theme = COLOR_SCHEMES.get(color_scheme, COLOR_SCHEMES["dark"])
    highlighted = selected_indices or ({selected_index} if selected_index is not None else set())
    projected = project_points(project, camera, width, height, mode=mode, center=center)
    along_edges = project.backbone_edges()
    _unused_along_edges, across_edges = project.all_edges()
    strand_colors = {strand.strand_id: strand.color for strand in project.strands}
    selected_next_along, selected_across = _selected_edge_sets(project, highlighted)
    effective_fit_zoom = fit_zoom or _fit_camera_zoom(project, camera, width, height, mode, center=center)
    effective_zoom = camera.zoom if camera.zoom > 0.0 else effective_fit_zoom
    radius = _scaled_marker_radius(sphere_size, effective_zoom, effective_fit_zoom)
    selected_radius = _scaled_selected_radius(radius)
    selected_backbone_width = backbone_width + 2
    svg = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        f'<rect x="0" y="0" width="{width}" height="{height}" fill="{theme["canvas_bg"]}"/>',
    ]
    for left, right in across_edges:
        x1, y1, _ = projected[left]
        x2, y2, _ = projected[right]
        dash_attr = _svg_dash_attr(_basepair_dash(project, left, right))
        svg.append(
            f'<line x1="{x1:.2f}" y1="{y1:.2f}" x2="{x2:.2f}" y2="{y2:.2f}" stroke="#64748b" stroke-width="{basepair_width}"{dash_attr}/>'
        )
    for left, right in along_edges:
        x1, y1, _ = projected[left]
        x2, y2, _ = projected[right]
        color = _backbone_color_for_base(project, left, strand_colors)
        svg.append(
            f'<line x1="{x1:.2f}" y1="{y1:.2f}" x2="{x2:.2f}" y2="{y2:.2f}" stroke="{color}" stroke-width="{backbone_width}"/>'
        )
    extension, _arrowshape, svg_arrow_length, svg_arrow_half_width = _strand_end_arrow_style(radius)
    for strand in project.strands:
        segment = _strand_end_arrow_segment(project, strand, projected, extension=extension)
        if segment is None:
            continue
        svg.extend(
            _svg_arrow_elements(
                segment,
                _strand_end_color(project, strand, strand_colors),
                backbone_width,
                arrow_length=svg_arrow_length,
                arrow_half_width=svg_arrow_half_width,
            )
        )
    if show_sticky_end:
        for group in project.sticky_end_groups():
            color = str(group["color"])
            for left_index, right_index in _sticky_end_pairs(group):
                x1, y1, _ = projected[left_index]
                x2, y2, _ = projected[right_index]
                svg.append(
                    f'<line x1="{x1:.2f}" y1="{y1:.2f}" x2="{x2:.2f}" y2="{y2:.2f}" stroke="{color}" stroke-width="1" stroke-dasharray="12,8"/>'
                )
            for indices in (group["left_indices"], group["right_indices"]):
                bounds = _projected_group_bounds(projected, indices, max(4.0, radius + 2.0))
                if bounds is None:
                    continue
                x0, y0, x1, y1 = bounds
                svg.append(
                    f'<rect x="{x0:.2f}" y="{y0:.2f}" width="{x1 - x0:.2f}" height="{y1 - y0:.2f}" fill="none" stroke="{color}" stroke-width="1"/>'
                )
    for left, right in selected_across:
        x1, y1, _ = projected[left]
        x2, y2, _ = projected[right]
        dash_attr = _svg_dash_attr(_basepair_dash(project, left, right))
        svg.append(
            f'<line x1="{x1:.2f}" y1="{y1:.2f}" x2="{x2:.2f}" y2="{y2:.2f}" stroke="{theme["selection_outline"]}" stroke-width="{basepair_width}"{dash_attr}/>'
        )
    for left, right in selected_next_along:
        x1, y1, _ = projected[left]
        x2, y2, _ = projected[right]
        strand_color = _backbone_color_for_base(project, left, strand_colors)
        svg.append(
            f'<line x1="{x1:.2f}" y1="{y1:.2f}" x2="{x2:.2f}" y2="{y2:.2f}" stroke="{theme["selection_outline"]}" stroke-width="{selected_backbone_width}"/>'
        )
        svg.append(
            f'<line x1="{x1:.2f}" y1="{y1:.2f}" x2="{x2:.2f}" y2="{y2:.2f}" stroke="{strand_color}" stroke-width="{backbone_width}"/>'
        )
    for index, (x, y, depth) in sorted(projected.items(), key=lambda item: item[1][2]):
        base = project.bases[index]
        is_selected = index in highlighted
        marker_radius = selected_radius if is_selected else radius
        outline = theme["selection_outline"] if is_selected else theme["base_outline"]
        width_attr = 2 if is_selected else 1
        fill = BASE_COLORS.get(base.nucleotide, BASE_COLORS[None])
        if str(getattr(base, "molecule", "DNA")).upper() == "RNA":
            svg.append(
                f'<rect x="{x - marker_radius:.2f}" y="{y - marker_radius:.2f}" width="{marker_radius * 2:.2f}" height="{marker_radius * 2:.2f}" '
                f'fill="{fill}" stroke="{outline}" stroke-width="{width_attr}"/>'
            )
        else:
            svg.append(
                f'<circle cx="{x:.2f}" cy="{y:.2f}" r="{marker_radius:.2f}" fill="{fill}" '
                f'stroke="{outline}" stroke-width="{width_attr}"/>'
            )
    svg.append("</svg>")
    return "\n".join(svg) + "\n"


def launch(project: TiamatProject | None = None) -> None:
    root = tk.Tk()
    TiamatViewer(root, project=project)
    root.mainloop()


def _project_vector(x: float, y: float, z: float, camera: Camera, mode: str) -> tuple[float, float, float]:
    if mode == "xy":
        return x, y, z
    if mode == "yz":
        return y, z, x
    if mode == "xz":
        return x, z, y
    return rotate_point(x, y, z, camera)


def _project_world_to_canvas(
    point: tuple[float, float, float],
    camera: Camera,
    width: int,
    height: int,
    mode: str,
    center: tuple[float, float, float],
    scale: float,
) -> tuple[float, float]:
    x, y, z = point
    rx, ry, _rz = _project_vector(x - center[0], y - center[1], z - center[2], camera, mode)
    return (
        width / 2.0 + camera.pan_x + rx * scale,
        height / 2.0 + camera.pan_y - ry * scale,
    )


def _project_base_marker_polygon(
    base: Base,
    camera: Camera,
    mode: str,
    center: tuple[float, float, float] | None,
    canvas_width: int,
    canvas_height: int,
    scale: float,
    radius: float,
) -> list[tuple[float, float]] | None:
    if base.plane_u is None or base.plane_v is None:
        return None
    world_radius = radius / max(scale, 1e-6)
    marker_points = _base_marker_world_points(base, world_radius)
    if marker_points is None:
        return None
    return [
        _project_world_to_canvas(point, camera, canvas_width, canvas_height, mode, center or (0.0, 0.0, 0.0), scale)
        for point in marker_points
    ]


def _base_marker_world_points(base: Base, world_radius: float) -> list[tuple[float, float, float]] | None:
    if base.plane_u is None or base.plane_v is None:
        return None
    plane_u = _normalize_vector(base.plane_u)
    plane_v = _normalize_vector(base.plane_v)
    if plane_u is None or plane_v is None:
        return None
    center = base.position
    if str(getattr(base, "molecule", "DNA")).upper() == "RNA":
        corners = (
            (-1.0, -1.0),
            (1.0, -1.0),
            (1.0, 1.0),
            (-1.0, 1.0),
        )
        return [
            _add_vectors(
                center,
                _add_vectors(
                    _scale_vector(plane_u, world_radius * corner_x),
                    _scale_vector(plane_v, world_radius * corner_y),
                ),
            )
            for corner_x, corner_y in corners
        ]
    points: list[tuple[float, float, float]] = []
    for step in range(16):
        angle = 2.0 * pi * step / 16.0
        points.append(
            _add_vectors(
                center,
                _add_vectors(
                    _scale_vector(plane_u, world_radius * cos(angle)),
                    _scale_vector(plane_v, world_radius * sin(angle)),
                ),
            )
        )
    return points


def _wrap_angle(value: float) -> float:
    return ((value + 180.0) % 360.0) - 180.0


def _copy_camera(camera: Camera) -> Camera:
    return Camera(
        rot_x=camera.rot_x,
        rot_y=camera.rot_y,
        rot_z=camera.rot_z,
        zoom=camera.zoom,
        pan_x=camera.pan_x,
        pan_y=camera.pan_y,
        orientation=_copy_orientation(camera.orientation),
    )


def _default_camera(mode: str) -> Camera:
    if mode == "3d":
        return Camera(rot_x=20.0, rot_y=-25.0, rot_z=0.0, zoom=0.0, pan_x=0.0, pan_y=0.0)
    return Camera(rot_x=0.0, rot_y=0.0, rot_z=0.0, zoom=0.0, pan_x=0.0, pan_y=0.0)


def _project_center(project: TiamatProject) -> tuple[float, float, float]:
    if not project.bases:
        return (0.0, 0.0, 0.0)
    xs = [base.x for base in project.bases.values()]
    ys = [base.y for base in project.bases.values()]
    zs = [base.z for base in project.bases.values()]
    return (
        (min(xs) + max(xs)) / 2.0,
        (min(ys) + max(ys)) / 2.0,
        (min(zs) + max(zs)) / 2.0,
    )


def _centered_coordinates(
    project: TiamatProject,
    center: tuple[float, float, float] | None = None,
) -> list[tuple[int, float, float, float]]:
    if not project.bases:
        return []
    origin = center or _project_center(project)
    return [
        (index, base.x - origin[0], base.y - origin[1], base.z - origin[2])
        for index, base in project.bases.items()
    ]


def _grid_axes_for_mode(mode: str) -> tuple[str, str]:
    if mode == "xy":
        return ("x", "y")
    if mode == "yz":
        return ("y", "z")
    return ("x", "z")


def _grid_world_point(mode: str, axis_x_value: float, axis_y_value: float) -> tuple[float, float, float]:
    if mode == "xy":
        return (axis_x_value, axis_y_value, 0.0)
    if mode == "yz":
        return (0.0, axis_x_value, axis_y_value)
    return (axis_x_value, 0.0, axis_y_value)


def _grid_fit_zoom(
    grid_extents: dict[str, float],
    width: int,
    height: int,
    mode: str,
) -> float:
    axis_x, axis_y = _grid_axes_for_mode(mode)
    span_x = max(grid_extents[f"{axis_x}_max"] - grid_extents[f"{axis_x}_min"], 1.0)
    span_y = max(grid_extents[f"{axis_y}_max"] - grid_extents[f"{axis_y}_min"], 1.0)
    padding = 48.0
    usable_width = max(width - padding, 1.0)
    usable_height = max(height - padding, 1.0)
    return min(usable_width / span_x, usable_height / span_y) * 0.92


def _grid_fit_zoom_for_mode(
    grid_extents: dict[str, float],
    camera: Camera,
    width: int,
    height: int,
    mode: str,
) -> float:
    if mode != "3d":
        return _grid_fit_zoom(grid_extents, width, height, mode)

    points = [
        (grid_extents["x_min"], grid_extents["y_min"], grid_extents["z_min"]),
        (grid_extents["x_min"], grid_extents["y_min"], grid_extents["z_max"]),
        (grid_extents["x_min"], grid_extents["y_max"], grid_extents["z_min"]),
        (grid_extents["x_min"], grid_extents["y_max"], grid_extents["z_max"]),
        (grid_extents["x_max"], grid_extents["y_min"], grid_extents["z_min"]),
        (grid_extents["x_max"], grid_extents["y_min"], grid_extents["z_max"]),
        (grid_extents["x_max"], grid_extents["y_max"], grid_extents["z_min"]),
        (grid_extents["x_max"], grid_extents["y_max"], grid_extents["z_max"]),
    ]
    projected = [_project_vector(x, y, z, camera, mode) for x, y, z in points]
    xs = [item[0] for item in projected]
    ys = [item[1] for item in projected]
    span_x = max(max(xs) - min(xs), 1.0)
    span_y = max(max(ys) - min(ys), 1.0)
    padding = 48.0
    usable_width = max(width - padding, 1.0)
    usable_height = max(height - padding, 1.0)
    return min(usable_width / span_x, usable_height / span_y) * 0.92


def _grid_values(min_value: float, max_value: float, spacing: float) -> list[float]:
    if spacing <= 0.0 or max_value < min_value:
        return []
    values: list[float] = []
    current = min_value
    tolerance = spacing * 1e-6
    while current <= max_value + tolerance:
        snapped = round(current / spacing) * spacing
        values.append(0.0 if abs(snapped) <= tolerance else round(snapped, 6))
        current += spacing
    return values


def _blend_hex(foreground: str, background: str, alpha: float) -> str:
    def parse(value: str) -> tuple[int, int, int]:
        value = value.lstrip("#")
        return (int(value[0:2], 16), int(value[2:4], 16), int(value[4:6], 16))

    fg_r, fg_g, fg_b = parse(foreground)
    bg_r, bg_g, bg_b = parse(background)
    alpha = max(0.0, min(alpha, 1.0))
    r = round(fg_r * alpha + bg_r * (1.0 - alpha))
    g = round(fg_g * alpha + bg_g * (1.0 - alpha))
    b = round(fg_b * alpha + bg_b * (1.0 - alpha))
    return f"#{r:02x}{g:02x}{b:02x}"


def _fit_camera_zoom(
    project: TiamatProject,
    camera: Camera,
    width: int,
    height: int,
    mode: str,
    center: tuple[float, float, float] | None = None,
) -> float:
    centered = _centered_coordinates(project, center=center)
    if not centered:
        return 1.0
    projected = [_project_vector(x, y, z, camera, mode) for _index, x, y, z in centered]
    xs = [item[0] for item in projected]
    ys = [item[1] for item in projected]
    span_x = max(max(xs) - min(xs), 1.0)
    span_y = max(max(ys) - min(ys), 1.0)
    padding = 48.0
    usable_width = max(width - padding, 1.0)
    usable_height = max(height - padding, 1.0)
    return min(usable_width / span_x, usable_height / span_y) * 0.92


def _scaled_marker_radius(base_size: float, zoom: float, fit_zoom: float) -> float:
    safe_fit_zoom = max(fit_zoom, 1e-6)
    safe_zoom = max(zoom, 1e-6)
    relative_scale = (safe_zoom / safe_fit_zoom) ** 0.5
    relative_scale = max(0.35, min(2.5, relative_scale))
    return max(1.2, min(48.0, base_size * relative_scale))


def _scaled_selected_radius(radius: float) -> float:
    return max(radius + 2.0, radius * 1.45)


def _copy_constraints_values(values: dict[str, dict[str, str]]) -> dict[str, dict[str, str]]:
    return {
        row_key: {column_key: str(value) for column_key, value in row_values.items()}
        for row_key, row_values in values.items()
    }


def _default_constraints_values() -> dict[str, dict[str, str]]:
    return _copy_constraints_values(DEFAULT_CONSTRAINT_VALUES)


def _blank_constraints_values() -> dict[str, dict[str, str]]:
    return {
        row_key: {column_key: "" for column_key, _column_label in CONSTRAINT_COLUMNS}
        for row_key, _row_label in CONSTRAINT_ROWS
    }


def _blank_custom_constraints_values() -> dict[str, dict[str, str]]:
    return {
        row_key: {column_key: "" for column_key, _column_label in CUSTOM_CONSTRAINT_COLUMNS}
        for row_key, _row_label in CONSTRAINT_ROWS
    }


def _normalize_constraints_values(
    values: dict[str, object] | None,
    columns: tuple[tuple[str, str], ...] = CONSTRAINT_COLUMNS,
    fallback: dict[str, dict[str, str]] | None = None,
) -> dict[str, dict[str, str]]:
    normalized = _copy_constraints_values(
        fallback
        or (_default_constraints_values() if columns == CONSTRAINT_COLUMNS else _blank_custom_constraints_values())
    )
    if not isinstance(values, dict):
        return normalized
    for row_key, _row_label in CONSTRAINT_ROWS:
        row_values = values.get(row_key)
        if not isinstance(row_values, dict):
            continue
        for column_key, _column_label in columns:
            raw_value = row_values.get(column_key, normalized[row_key][column_key])
            normalized[row_key][column_key] = "" if raw_value is None else str(raw_value).strip()
    _synchronize_helicity_values(normalized, columns, strict=False)
    return normalized


def _validate_constraints_values(
    values: dict[str, dict[str, str]],
    columns: tuple[tuple[str, str], ...] = CONSTRAINT_COLUMNS,
) -> dict[str, dict[str, str]]:
    normalized = _normalize_constraints_values(
        values,
        columns=columns,
        fallback=_blank_constraints_values() if columns == CONSTRAINT_COLUMNS else _blank_custom_constraints_values(),
    )
    _synchronize_helicity_values(normalized, columns, strict=True)
    for row_key, row_label in CONSTRAINT_ROWS:
        for column_key, column_label in columns:
            text = normalized[row_key][column_key]
            if text == "":
                continue
            try:
                float(text)
            except ValueError as exc:
                raise ValueError(f"{row_label} / {column_label} must be numeric or blank.") from exc
    return normalized


def _project_constraints(project: TiamatProject | None) -> dict[str, dict[str, str]]:
    if project is None:
        return _default_constraints_values()
    return _normalize_constraints_values(project.metadata.get("constraints"))


def _project_custom_constraints(project: TiamatProject | None) -> dict[str, dict[str, str]]:
    if project is None:
        return _blank_custom_constraints_values()
    return _normalize_constraints_values(
        project.metadata.get("custom_constraints"),
        columns=CUSTOM_CONSTRAINT_COLUMNS,
        fallback=_blank_custom_constraints_values(),
    )


def _format_rotation_value(value: float) -> str:
    return f"{value:.4f}".rstrip("0").rstrip(".")


def _format_helicity_value(value: float) -> str:
    return f"{value:.2f}"


def _helicity_from_rotation_value(rotation: float) -> float:
    if abs(rotation) <= 1e-9:
        raise ValueError("Rotation per Base Pair must be non-zero.")
    return -360.0 / rotation


def _rotation_from_helicity_value(helicity: float) -> float:
    if abs(helicity) <= 1e-9:
        raise ValueError("Helicity must be non-zero.")
    return -360.0 / helicity


def _synchronize_helicity_values(
    values: dict[str, dict[str, str]],
    columns: tuple[tuple[str, str], ...],
    strict: bool,
) -> None:
    if "rotation_per_base_pair_deg" not in values or "helicity_bp_per_turn" not in values:
        return
    for column_key, _column_label in columns:
        rotation_text = values["rotation_per_base_pair_deg"].get(column_key, "").strip()
        helicity_text = values["helicity_bp_per_turn"].get(column_key, "").strip()
        if rotation_text:
            try:
                rotation = float(rotation_text)
                values["helicity_bp_per_turn"][column_key] = _format_helicity_value(_helicity_from_rotation_value(rotation))
            except ValueError:
                if strict:
                    raise
        elif helicity_text:
            try:
                helicity = float(helicity_text)
                values["rotation_per_base_pair_deg"][column_key] = _format_rotation_value(_rotation_from_helicity_value(helicity))
                values["helicity_bp_per_turn"][column_key] = _format_helicity_value(helicity)
            except ValueError:
                if strict:
                    raise


def _view_scale(
    project: TiamatProject | None,
    grid_extents: dict[str, float],
    state: ViewportState,
    width: int,
    height: int,
    center: tuple[float, float, float],
) -> float:
    if state.camera.zoom > 0.0:
        return state.camera.zoom
    if project and project.bases:
        return _fit_camera_zoom(project, state.camera, width, height, state.mode, center=center)
    return _grid_fit_zoom_for_mode(grid_extents, state.camera, width, height, state.mode)


def _canvas_to_world(
    x: float,
    y: float,
    camera: Camera,
    width: int,
    height: int,
    mode: str,
    center: tuple[float, float, float],
    scale: float,
) -> tuple[float, float, float]:
    safe_scale = max(scale, 1e-6)
    horizontal = (x - width / 2.0 - camera.pan_x) / safe_scale
    vertical = -(y - height / 2.0 - camera.pan_y) / safe_scale
    if mode == "xy":
        return (center[0] + horizontal, center[1] + vertical, center[2])
    if mode == "yz":
        return (center[0], center[1] + horizontal, center[2] + vertical)
    if mode == "xz":
        return (center[0] + horizontal, center[1], center[2] + vertical)
    right_axis, up_axis, _forward_axis = _camera_world_axes(camera)
    return (
        center[0] + right_axis[0] * horizontal + up_axis[0] * vertical,
        center[1] + right_axis[1] * horizontal + up_axis[1] * vertical,
        center[2] + right_axis[2] * horizontal + up_axis[2] * vertical,
    )


def _snap_create_helix_endpoint(
    start_canvas: tuple[int, int],
    end_canvas: tuple[int, int],
    camera: Camera,
    width: int,
    height: int,
    mode: str,
    center: tuple[float, float, float],
    scale: float,
) -> tuple[int, int]:
    start_world = _canvas_to_world(start_canvas[0], start_canvas[1], camera, width, height, mode, center, scale)
    end_world = _canvas_to_world(end_canvas[0], end_canvas[1], camera, width, height, mode, center, scale)
    vector = (
        end_world[0] - start_world[0],
        end_world[1] - start_world[1],
        end_world[2] - start_world[2],
    )
    if max(abs(vector[0]), abs(vector[1]), abs(vector[2])) <= 1e-9:
        return end_canvas
    axis_index = max(range(3), key=lambda index: abs(vector[index]))
    axis_vector = ((1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0))[axis_index]
    signed_length = vector[axis_index]
    snapped_world = (
        start_world[0] + axis_vector[0] * signed_length,
        start_world[1] + axis_vector[1] * signed_length,
        start_world[2] + axis_vector[2] * signed_length,
    )
    snapped_canvas = _project_world_to_canvas(snapped_world, camera, width, height, mode, center, scale)
    return (int(round(snapped_canvas[0])), int(round(snapped_canvas[1])))


def _estimate_created_base_count(axis_length: float, constraints_values: dict[str, dict[str, str]]) -> int:
    try:
        rise = _constraint_value(constraints_values, "rise_per_base_pair_nm", "b_median")
    except ValueError:
        return DEFAULT_CREATE_BASE_COUNT
    if axis_length <= 1e-6 or rise <= 1e-6:
        return DEFAULT_CREATE_BASE_COUNT
    return max(1, round(axis_length / rise))


def _custom_constraints_defined(constraints_values: dict[str, dict[str, str]]) -> bool:
    try:
        for row_key, _row_label in CONSTRAINT_ROWS:
            _constraint_value(constraints_values, row_key, "x_median")
    except ValueError:
        return False
    return True


def _create_helix_molecule_option(molecule: str, family: str) -> str:
    normalized_molecule = molecule.strip().upper()
    family_mapping = {"A-TYPE": "A", "B-TYPE": "B", "X-TYPE": "X"}
    normalized_family = family_mapping.get(family.strip().upper())
    if normalized_family is None:
        raise ValueError(f"Unsupported helix structure: {family}.")
    return f"{normalized_molecule} {normalized_family}"


def _constraint_value(
    constraints_values: dict[str, dict[str, str]],
    row_key: str,
    column_key: str,
) -> float:
    try:
        row_values = constraints_values[row_key]
    except KeyError as exc:
        raise ValueError(f"Missing constraint row: {row_key}.") from exc
    raw_value = str(row_values.get(column_key, "")).strip()
    if raw_value == "":
        raise ValueError(f"Constraint {row_key} / {column_key} is blank.")
    try:
        return float(raw_value)
    except ValueError as exc:
        raise ValueError(f"Constraint {row_key} / {column_key} must be numeric.") from exc


def _constraint_profile_for_molecule_option(
    molecule_option: str,
    constraints_values: dict[str, dict[str, str]],
) -> dict[str, object]:
    option = molecule_option.strip().upper()
    mapping = {
        "DNA-DNA B": (("DNA", "DNA"), "b"),
        "DNA-DNA A": (("DNA", "DNA"), "a"),
        "DNA-DNA X": (("DNA", "DNA"), "x"),
        "DNA-RNA B": (("DNA", "RNA"), "b"),
        "DNA-RNA A": (("DNA", "RNA"), "a"),
        "DNA-RNA X": (("DNA", "RNA"), "x"),
        "RNA-DNA B": (("RNA", "DNA"), "b"),
        "RNA-DNA A": (("RNA", "DNA"), "a"),
        "RNA-DNA X": (("RNA", "DNA"), "x"),
        "RNA-RNA B": (("RNA", "RNA"), "b"),
        "RNA-RNA A": (("RNA", "RNA"), "a"),
        "RNA-RNA X": (("RNA", "RNA"), "x"),
    }
    if option not in mapping:
        raise ValueError(f"Unsupported molecule option: {molecule_option}.")
    strand_molecules, helix_family = mapping[option]
    median_column = {"b": "b_median", "a": "a_median", "x": "x_median"}[helix_family]
    return {
        "strand_molecules": strand_molecules,
        "helix_family": helix_family,
        "rise": _constraint_value(constraints_values, "rise_per_base_pair_nm", median_column),
        "rotation": _constraint_value(constraints_values, "rotation_per_base_pair_deg", median_column),
        "diameter": _constraint_value(constraints_values, "diameter_nm", median_column),
        "inclination": _constraint_value(constraints_values, "inclination_deg", median_column),
        "minor_groove_angle": _constraint_value(constraints_values, "minor_groove_angle_deg", median_column),
    }


def _append_created_structure(
    project: TiamatProject,
    start_point: tuple[float, float, float],
    direction: tuple[float, float, float],
    count: int,
    backbone_rotation: float,
    molecule_option: str,
    structure_type: str,
    constraints_values: dict[str, dict[str, str]],
) -> list[int]:
    created = _build_created_structure(
        start_point=start_point,
        direction=direction,
        count=count,
        backbone_rotation=backbone_rotation,
        molecule_option=molecule_option,
        structure_type=structure_type,
        constraints_values=constraints_values,
        start_index=max(project.bases, default=0) + 1,
    )
    project.bases.update(created)
    project.refresh_strands()
    return sorted(created)


def _base_terminal_roles(base: Base) -> frozenset[str]:
    roles: set[str] = set()
    if base.up is None:
        roles.add("5")
    if base.down is None:
        roles.add("3")
    return frozenset(roles)


def _estimate_free_strand_base_count(
    start_point: tuple[float, float, float],
    end_point: tuple[float, float, float],
    spacing_nm: float = 0.6,
) -> int:
    distance = sqrt(
        (end_point[0] - start_point[0]) ** 2
        + (end_point[1] - start_point[1]) ** 2
        + (end_point[2] - start_point[2]) ** 2
    )
    if spacing_nm <= 1e-9:
        return 1
    return max(1, round(distance / spacing_nm))


def _resolve_free_strand_creation(
    first: FreeStrandEndpoint,
    second: FreeStrandEndpoint,
) -> dict[str, object]:
    if first.base_index is not None and second.base_index is not None and first.base_index == second.base_index:
        raise ValueError("Select a different second endpoint.")

    if first.base_index is None and second.base_index is None:
        return {
            "start_point": first.point,
            "end_point": second.point,
            "start_anchor": None,
            "end_anchor": None,
            "molecule": "DNA",
        }

    if first.base_index is not None and second.base_index is not None:
        possible = []
        if "3" in first.terminal_roles and "5" in second.terminal_roles:
            possible.append((first, second))
        if "3" in second.terminal_roles and "5" in first.terminal_roles:
            possible.append((second, first))
        if not possible:
            raise ValueError("Free strand must connect 3' and 5' ends.")
        start_anchor, end_anchor = possible[0]
        return {
            "start_point": start_anchor.point,
            "end_point": end_anchor.point,
            "start_anchor": start_anchor.base_index,
            "end_anchor": end_anchor.base_index,
            "molecule": str(start_anchor.molecule or end_anchor.molecule or "DNA").upper(),
        }

    terminal = first if first.base_index is not None else second
    empty = second if terminal is first else first
    if terminal.terminal_roles == frozenset({"3"}):
        return {
            "start_point": terminal.point,
            "end_point": empty.point,
            "start_anchor": terminal.base_index,
            "end_anchor": None,
            "molecule": str(terminal.molecule or "DNA").upper(),
        }
    if terminal.terminal_roles == frozenset({"5"}):
        return {
            "start_point": empty.point,
            "end_point": terminal.point,
            "start_anchor": None,
            "end_anchor": terminal.base_index,
            "molecule": str(terminal.molecule or "DNA").upper(),
        }
    raise ValueError("Ambiguous terminal end. Free strand creation requires a distinct 3' or 5' end.")


def _append_free_strand(
    project: TiamatProject,
    start_point: tuple[float, float, float],
    end_point: tuple[float, float, float],
    count: int,
    molecule: str,
    start_anchor: int | None = None,
    end_anchor: int | None = None,
) -> list[int]:
    if count <= 0:
        raise ValueError("Free strand length must be greater than 0.")
    if start_anchor is not None:
        start_base = project.bases[start_anchor]
        if start_base.down is not None:
            raise ValueError("Selected start base is not a 3' terminal base.")
    if end_anchor is not None:
        end_base = project.bases[end_anchor]
        if end_base.up is not None:
            raise ValueError("Selected end base is not a 5' terminal base.")

    delta = (
        end_point[0] - start_point[0],
        end_point[1] - start_point[1],
        end_point[2] - start_point[2],
    )
    direction = _normalize_vector(delta) or (1.0, 0.0, 0.0)
    plane_u, plane_v = _perpendicular_basis(direction)
    start_inclusive = start_anchor is None
    end_inclusive = end_anchor is None
    if start_inclusive and end_inclusive:
        denominator = count - 1
        offset_start = 0
    elif start_inclusive or end_inclusive:
        denominator = count
        offset_start = 0 if start_inclusive else 1
    else:
        denominator = count + 1
        offset_start = 1

    start_index = max(project.bases, default=0) + 1
    created_indices: list[int] = []
    previous_index = start_anchor
    for offset in range(count):
        if denominator <= 0:
            fraction = 0.5
        else:
            fraction = (offset_start + offset) / denominator
        x = start_point[0] + delta[0] * fraction
        y = start_point[1] + delta[1] * fraction
        z = start_point[2] + delta[2] * fraction
        index = start_index + offset
        created_indices.append(index)
        project.bases[index] = Base(
            index=index,
            object_id=str(index),
            x=x,
            y=y,
            z=z,
            molecule=str(molecule).upper(),
            up=previous_index,
            nucleotide=None,
            plane_u=plane_u,
            plane_v=plane_v,
        )
        if previous_index is not None:
            project.bases[previous_index].down = index
        previous_index = index

    if end_anchor is not None:
        project.bases[previous_index].down = end_anchor
        project.bases[end_anchor].up = previous_index

    project.refresh_strands()
    return created_indices


def _build_created_structure(
    start_point: tuple[float, float, float],
    direction: tuple[float, float, float],
    count: int,
    backbone_rotation: float,
    molecule_option: str,
    structure_type: str,
    constraints_values: dict[str, dict[str, str]],
    start_index: int,
) -> dict[int, Base]:
    if count <= 0:
        raise ValueError("Number of Bases must be greater than 0.")
    axis = _normalize_vector(direction)
    if axis is None:
        axis = (1.0, 0.0, 0.0)
    profile = _constraint_profile_for_molecule_option(molecule_option, constraints_values)
    first_molecule, second_molecule = profile["strand_molecules"]
    rise = float(profile["rise"])
    rotation = float(profile["rotation"])
    radius = max(float(profile["diameter"]) / 2.0, 0.0)
    inclination = float(profile["inclination"])
    minor_groove_angle = float(profile["minor_groove_angle"])
    structure_kind = structure_type.strip().lower()
    if structure_kind not in {"helix", "strand"}:
        raise ValueError(f"Unsupported structure type: {structure_type}.")

    radial_u, radial_v = _perpendicular_basis(axis)
    inclination_offset = radius * tan(radians(inclination))
    created: dict[int, Base] = {}
    first_indices: list[int] = []
    second_indices: list[int] = []
    next_index = start_index

    for position in range(count):
        axis_offset = _scale_vector(axis, rise * position)
        axis_point = _add_vectors(start_point, axis_offset)
        phase = radians(backbone_rotation - rotation * position)
        first_radial = _radial_vector(radial_u, radial_v, phase)
        second_radial = _radial_vector(radial_u, radial_v, phase + radians(minor_groove_angle))
        first_position = _add_vectors(
            _add_vectors(axis_point, _scale_vector(first_radial, radius)),
            _scale_vector(axis, inclination_offset),
        )
        pair_frame = _base_plane_frame(
            axis=axis,
            first_radial=first_radial,
            second_radial=second_radial,
            inclination_deg=inclination,
        )
        first_index = next_index
        next_index += 1
        first_indices.append(first_index)
        created[first_index] = Base(
            index=first_index,
            object_id=str(first_index),
            x=first_position[0],
            y=first_position[1],
            z=first_position[2],
            molecule=first_molecule,
            nucleotide=None,
            plane_u=pair_frame[0],
            plane_v=pair_frame[1],
        )
        if structure_kind == "helix":
            second_position = _add_vectors(
                _add_vectors(axis_point, _scale_vector(second_radial, radius)),
                _scale_vector(axis, -inclination_offset),
            )
            second_index = next_index
            next_index += 1
            second_indices.append(second_index)
            created[second_index] = Base(
                index=second_index,
                object_id=str(second_index),
                x=second_position[0],
                y=second_position[1],
                z=second_position[2],
                molecule=second_molecule,
                nucleotide=None,
                plane_u=pair_frame[0],
                plane_v=pair_frame[1],
            )

    for position, index in enumerate(first_indices):
        created[index].up = first_indices[position - 1] if position > 0 else None
        created[index].down = first_indices[position + 1] if position < len(first_indices) - 1 else None
        if structure_kind == "helix":
            created[index].across = second_indices[position]

    if structure_kind == "helix":
        for position, index in enumerate(second_indices):
            created[index].up = second_indices[position + 1] if position < len(second_indices) - 1 else None
            created[index].down = second_indices[position - 1] if position > 0 else None
            created[index].across = first_indices[position]

    return created


def _normalize_vector(vector: tuple[float, float, float]) -> tuple[float, float, float] | None:
    magnitude = sqrt(vector[0] * vector[0] + vector[1] * vector[1] + vector[2] * vector[2])
    if magnitude <= 1e-9:
        return None
    return (vector[0] / magnitude, vector[1] / magnitude, vector[2] / magnitude)


def _add_vectors(left: tuple[float, float, float], right: tuple[float, float, float]) -> tuple[float, float, float]:
    return (left[0] + right[0], left[1] + right[1], left[2] + right[2])


def _scale_vector(vector: tuple[float, float, float], scale: float) -> tuple[float, float, float]:
    return (vector[0] * scale, vector[1] * scale, vector[2] * scale)


def _subtract_vectors(left: tuple[float, float, float], right: tuple[float, float, float]) -> tuple[float, float, float]:
    return (left[0] - right[0], left[1] - right[1], left[2] - right[2])


def _dot_product(left: tuple[float, float, float], right: tuple[float, float, float]) -> float:
    return left[0] * right[0] + left[1] * right[1] + left[2] * right[2]


def _cross_product(left: tuple[float, float, float], right: tuple[float, float, float]) -> tuple[float, float, float]:
    return (
        left[1] * right[2] - left[2] * right[1],
        left[2] * right[0] - left[0] * right[2],
        left[0] * right[1] - left[1] * right[0],
    )


def _perpendicular_basis(
    axis: tuple[float, float, float],
) -> tuple[tuple[float, float, float], tuple[float, float, float]]:
    reference = (0.0, 0.0, 1.0) if abs(axis[2]) < 0.9 else (0.0, 1.0, 0.0)
    first = _normalize_vector(_cross_product(axis, reference))
    if first is None:
        first = (1.0, 0.0, 0.0)
    second = _normalize_vector(_cross_product(axis, first))
    if second is None:
        second = (0.0, 1.0, 0.0)
    return first, second


def _radial_vector(
    radial_u: tuple[float, float, float],
    radial_v: tuple[float, float, float],
    angle_radians: float,
) -> tuple[float, float, float]:
    return _add_vectors(
        _scale_vector(radial_u, cos(angle_radians)),
        _scale_vector(radial_v, sin(angle_radians)),
    )


def _base_plane_frame(
    axis: tuple[float, float, float],
    first_radial: tuple[float, float, float],
    second_radial: tuple[float, float, float],
    inclination_deg: float,
) -> tuple[tuple[float, float, float], tuple[float, float, float], tuple[float, float, float]]:
    axis_unit = _normalize_vector(axis) or (1.0, 0.0, 0.0)
    midpoint_direction = _normalize_vector(_add_vectors(first_radial, second_radial)) or first_radial
    inclination_radians = radians(inclination_deg)
    plane_normal = _normalize_vector(
        _add_vectors(
            _scale_vector(axis_unit, cos(inclination_radians)),
            _scale_vector(midpoint_direction, sin(inclination_radians)),
        )
    ) or axis_unit
    pair_direction = _normalize_vector(_subtract_vectors(second_radial, first_radial)) or midpoint_direction
    in_plane = _subtract_vectors(pair_direction, _scale_vector(plane_normal, _dot_product(pair_direction, plane_normal)))
    plane_u = _normalize_vector(in_plane)
    if plane_u is None:
        plane_u, _unused_plane_v = _perpendicular_basis(plane_normal)
    plane_v = _normalize_vector(_cross_product(plane_normal, plane_u))
    if plane_v is None:
        plane_v = _normalize_vector(_cross_product(plane_normal, midpoint_direction)) or (0.0, 1.0, 0.0)
    return plane_u, plane_v, plane_normal


def _clipboard_render_scale(width: int, height: int, max_dimension: float = CLIPBOARD_RENDER_MAX_DIMENSION) -> float:
    longest_edge = max(float(width), float(height), 1.0)
    target = max(float(max_dimension), 1.0)
    return max(1.0, target / longest_edge)


def _copy_orientation(
    orientation: tuple[tuple[float, float, float], tuple[float, float, float], tuple[float, float, float]] | None,
) -> tuple[tuple[float, float, float], tuple[float, float, float], tuple[float, float, float]] | None:
    if orientation is None:
        return None
    return tuple(tuple(value for value in row) for row in orientation)


def _apply_orientation(
    orientation: tuple[tuple[float, float, float], tuple[float, float, float], tuple[float, float, float]],
    x: float,
    y: float,
    z: float,
) -> tuple[float, float, float]:
    return (
        orientation[0][0] * x + orientation[0][1] * y + orientation[0][2] * z,
        orientation[1][0] * x + orientation[1][1] * y + orientation[1][2] * z,
        orientation[2][0] * x + orientation[2][1] * y + orientation[2][2] * z,
    )


def _rotation_matrix_x(angle_degrees: float) -> tuple[tuple[float, float, float], tuple[float, float, float], tuple[float, float, float]]:
    angle = radians(angle_degrees)
    c = cos(angle)
    s = sin(angle)
    return (
        (1.0, 0.0, 0.0),
        (0.0, c, -s),
        (0.0, s, c),
    )


def _rotation_matrix_y(angle_degrees: float) -> tuple[tuple[float, float, float], tuple[float, float, float], tuple[float, float, float]]:
    angle = radians(angle_degrees)
    c = cos(angle)
    s = sin(angle)
    return (
        (c, 0.0, s),
        (0.0, 1.0, 0.0),
        (-s, 0.0, c),
    )


def _rotation_matrix_z(angle_degrees: float) -> tuple[tuple[float, float, float], tuple[float, float, float], tuple[float, float, float]]:
    angle = radians(angle_degrees)
    c = cos(angle)
    s = sin(angle)
    return (
        (c, -s, 0.0),
        (s, c, 0.0),
        (0.0, 0.0, 1.0),
    )


def _matrix_multiply(
    left: tuple[tuple[float, float, float], tuple[float, float, float], tuple[float, float, float]],
    right: tuple[tuple[float, float, float], tuple[float, float, float], tuple[float, float, float]],
) -> tuple[tuple[float, float, float], tuple[float, float, float], tuple[float, float, float]]:
    result: list[tuple[float, float, float]] = []
    for row in range(3):
        values = []
        for column in range(3):
            values.append(
                left[row][0] * right[0][column]
                + left[row][1] * right[1][column]
                + left[row][2] * right[2][column]
            )
        result.append((values[0], values[1], values[2]))
    return (result[0], result[1], result[2])


def _camera_orientation(
    camera: Camera,
) -> tuple[tuple[float, float, float], tuple[float, float, float], tuple[float, float, float]]:
    if camera.orientation is not None:
        return camera.orientation
    return _matrix_multiply(
        _rotation_matrix_z(camera.rot_z),
        _matrix_multiply(_rotation_matrix_y(camera.rot_y), _rotation_matrix_x(camera.rot_x)),
    )


def _apply_trackball_rotation(
    camera: Camera,
    dx: float,
    dy: float,
    sensitivity: float = 0.45,
) -> tuple[tuple[float, float, float], tuple[float, float, float], tuple[float, float, float]]:
    base_orientation = _camera_orientation(camera)
    delta_rotation = _matrix_multiply(
        _rotation_matrix_y(dx * sensitivity),
        _rotation_matrix_x(-dy * sensitivity),
    )
    return _matrix_multiply(delta_rotation, base_orientation)


def _selection_status_text(
    project: TiamatProject,
    selected_indices: set[int],
    primary_selected: int | None,
) -> str:
    if not selected_indices:
        return "Selected: -"
    if len(selected_indices) == 1:
        if primary_selected in selected_indices:
            return f"Selected base ID: {primary_selected}"
        return f"Selected base ID: {min(selected_indices)}"

    summary = f"Selected bases: {len(selected_indices)}"
    if len(selected_indices) == 2:
        first, second = sorted(selected_indices)
        summary += f" | Distance: {_distance_between_bases_nm(project, first, second):.2f} nm"
    return summary


def _distance_between_bases_nm(project: TiamatProject, left: int, right: int) -> float:
    first = project.bases[left]
    second = project.bases[right]
    dx = first.x - second.x
    dy = first.y - second.y
    dz = first.z - second.z
    return (dx * dx + dy * dy + dz * dz) ** 0.5


def _selection_center(project: TiamatProject, selected_indices: set[int]) -> tuple[float, float, float]:
    count = max(len(selected_indices), 1)
    total_x = 0.0
    total_y = 0.0
    total_z = 0.0
    for index in selected_indices:
        base = project.bases[index]
        total_x += base.x
        total_y += base.y
        total_z += base.z
    return (total_x / count, total_y / count, total_z / count)


def _projected_selection_center(
    projected: dict[int, tuple[float, float, float]],
    selected_indices: set[int],
) -> tuple[float, float] | None:
    points = [projected[index] for index in selected_indices if index in projected]
    if not points:
        return None
    x = sum(point[0] for point in points) / len(points)
    y = sum(point[1] for point in points) / len(points)
    return (x, y)


def _apply_translation_to_project(
    project: TiamatProject,
    original_positions: dict[int, tuple[float, float, float]],
    translation: tuple[float, float, float],
) -> None:
    dx, dy, dz = translation
    for index, (x, y, z) in original_positions.items():
        base = project.bases[index]
        base.x = x + dx
        base.y = y + dy
        base.z = z + dz


def _apply_rotation_to_project(
    project: TiamatProject,
    original_positions: dict[int, tuple[float, float, float]],
    center: tuple[float, float, float],
    rotation: tuple[tuple[float, float, float], tuple[float, float, float], tuple[float, float, float]],
) -> None:
    center_x, center_y, center_z = center
    for index, (x, y, z) in original_positions.items():
        rotated_x, rotated_y, rotated_z = _apply_orientation(
            rotation,
            x - center_x,
            y - center_y,
            z - center_z,
        )
        base = project.bases[index]
        base.x = center_x + rotated_x
        base.y = center_y + rotated_y
        base.z = center_z + rotated_z
        if base.plane_u is not None:
            base.plane_u = _apply_orientation(rotation, *base.plane_u)
        if base.plane_v is not None:
            base.plane_v = _apply_orientation(rotation, *base.plane_v)


def _drag_translation_vector(
    mode: str,
    camera: Camera,
    dx: float,
    dy: float,
    scale: float,
    constrain_axis: bool = False,
) -> tuple[float, float, float]:
    scale = max(scale, 1e-6)
    horizontal = dx / scale
    vertical = -dy / scale
    if constrain_axis:
        if abs(horizontal) >= abs(vertical):
            vertical = 0.0
        else:
            horizontal = 0.0
    if mode == "xy":
        return (horizontal, vertical, 0.0)
    if mode == "yz":
        return (0.0, horizontal, vertical)
    if mode == "xz":
        return (horizontal, 0.0, vertical)
    right_axis, up_axis, _forward_axis = _camera_world_axes(camera)
    return (
        right_axis[0] * horizontal + up_axis[0] * vertical,
        right_axis[1] * horizontal + up_axis[1] * vertical,
        right_axis[2] * horizontal + up_axis[2] * vertical,
    )


def _camera_world_axes(camera: Camera) -> tuple[tuple[float, float, float], tuple[float, float, float], tuple[float, float, float]]:
    orientation = _camera_orientation(camera)
    right_axis = (orientation[0][0], orientation[1][0], orientation[2][0])
    up_axis = (orientation[0][1], orientation[1][1], orientation[2][1])
    forward_axis = (orientation[0][2], orientation[1][2], orientation[2][2])
    return right_axis, up_axis, forward_axis


def _drag_rotation_matrix_planar(
    mode: str,
    dx: float,
    dy: float,
    sensitivity: float = 0.45,
    constrain_axis: bool = False,
) -> tuple[tuple[float, float, float], tuple[float, float, float], tuple[float, float, float]]:
    axis_names = {
        "xy": ("y", "x"),
        "yz": ("z", "y"),
        "xz": ("z", "x"),
    }.get(mode)
    if axis_names is None:
        return (
            (1.0, 0.0, 0.0),
            (0.0, 1.0, 0.0),
            (0.0, 0.0, 1.0),
        )
    horizontal_drag_axis = _world_axis(axis_names[0])
    vertical_drag_axis = _world_axis(axis_names[1])
    horizontal_rotation = _rotation_matrix_axis(horizontal_drag_axis, dx * sensitivity)
    vertical_rotation = _rotation_matrix_axis(vertical_drag_axis, -dy * sensitivity)
    if constrain_axis:
        return horizontal_rotation if abs(dx) >= abs(dy) else vertical_rotation
    return _matrix_multiply(horizontal_rotation, vertical_rotation)


def _drag_rotation_matrix_3d(
    camera: Camera,
    dx: float,
    dy: float,
    sensitivity: float = 0.45,
    constrain_axis: bool = False,
) -> tuple[tuple[float, float, float], tuple[float, float, float], tuple[float, float, float]]:
    right_axis, up_axis, _forward_axis = _camera_world_axes(camera)
    horizontal = _rotation_matrix_axis(up_axis, dx * sensitivity)
    vertical = _rotation_matrix_axis(right_axis, -dy * sensitivity)
    if constrain_axis:
        return horizontal if abs(dx) >= abs(dy) else vertical
    return _matrix_multiply(horizontal, vertical)


def _rotation_matrix_axis(
    axis: tuple[float, float, float],
    angle_degrees: float,
) -> tuple[tuple[float, float, float], tuple[float, float, float], tuple[float, float, float]]:
    axis_x, axis_y, axis_z = axis
    magnitude = sqrt(axis_x * axis_x + axis_y * axis_y + axis_z * axis_z)
    if magnitude <= 1e-9 or abs(angle_degrees) <= 1e-9:
        return (
            (1.0, 0.0, 0.0),
            (0.0, 1.0, 0.0),
            (0.0, 0.0, 1.0),
        )
    axis_x /= magnitude
    axis_y /= magnitude
    axis_z /= magnitude
    angle = radians(angle_degrees)
    c = cos(angle)
    s = sin(angle)
    t = 1.0 - c
    return (
        (t * axis_x * axis_x + c, t * axis_x * axis_y - s * axis_z, t * axis_x * axis_z + s * axis_y),
        (t * axis_x * axis_y + s * axis_z, t * axis_y * axis_y + c, t * axis_y * axis_z - s * axis_x),
        (t * axis_x * axis_z - s * axis_y, t * axis_y * axis_z + s * axis_x, t * axis_z * axis_z + c),
    )


def _world_axis(name: str) -> tuple[float, float, float]:
    if name == "x":
        return (1.0, 0.0, 0.0)
    if name == "y":
        return (0.0, 1.0, 0.0)
    return (0.0, 0.0, 1.0)


def _strand_end_arrow_segment(
    project: TiamatProject,
    strand,
    projected: dict[int, tuple[float, float, float]],
    extension: float = 18.0,
) -> tuple[float, float, float, float] | None:
    if len(strand.base_indices) < 2:
        return None
    end_index = strand.base_indices[-1]
    if project.bases[end_index].down is not None:
        return None
    previous_index = strand.base_indices[-2]
    if previous_index not in projected or end_index not in projected:
        return None
    x1, y1, _ = projected[previous_index]
    x2, y2, _ = projected[end_index]
    dx = x2 - x1
    dy = y2 - y1
    length = (dx * dx + dy * dy) ** 0.5
    if length <= 1e-6:
        return None
    ux = dx / length
    uy = dy / length
    return (x2, y2, x2 + ux * extension, y2 + uy * extension)


def _strand_end_arrow_style(radius: float) -> tuple[float, tuple[float, float, float], float, float]:
    extension = max(6.0, radius * 4.5)
    scale = extension / 18.0
    arrowshape = (
        max(3.0, 7.0 * scale),
        max(4.0, 9.0 * scale),
        max(1.5, 3.0 * scale),
    )
    arrow_length = max(3.5, 7.0 * scale)
    arrow_half_width = max(1.5, 3.0 * scale)
    return extension, arrowshape, arrow_length, arrow_half_width


def _svg_arrow_elements(
    segment: tuple[float, float, float, float],
    color: str,
    width: int,
    arrow_length: float = 7.0,
    arrow_half_width: float = 3.0,
) -> list[str]:
    x1, y1, x2, y2 = segment
    dx = x2 - x1
    dy = y2 - y1
    length = (dx * dx + dy * dy) ** 0.5
    if length <= 1e-6:
        return []
    ux = dx / length
    uy = dy / length
    base_x = x2 - ux * arrow_length
    base_y = y2 - uy * arrow_length
    perp_x = -uy * arrow_half_width
    perp_y = ux * arrow_half_width
    return [
        f'<line x1="{x1:.2f}" y1="{y1:.2f}" x2="{base_x:.2f}" y2="{base_y:.2f}" stroke="{color}" stroke-width="{width}"/>',
        f'<polygon points="{x2:.2f},{y2:.2f} {base_x + perp_x:.2f},{base_y + perp_y:.2f} {base_x - perp_x:.2f},{base_y - perp_y:.2f}" fill="{color}"/>',
    ]


def _backbone_color_for_base(
    project: TiamatProject,
    start_index: int,
    strand_colors: dict[int, str],
) -> str:
    base = project.bases[start_index]
    if getattr(base, "backbone_color", None):
        return str(base.backbone_color)
    strand_id = base.strand_id
    return strand_colors.get(strand_id, "#cbd5e1")


def _sticky_end_pairs(group: dict[str, object]) -> list[tuple[int, int]]:
    left_indices = [int(value) for value in group.get("left_indices", [])]
    right_indices = [int(value) for value in group.get("right_indices", [])]
    return list(zip(left_indices, reversed(right_indices)))


def _projected_group_bounds(
    projected: dict[int, tuple[float, float, float]],
    indices: list[int],
    padding: float,
) -> tuple[float, float, float, float] | None:
    points = [projected[index] for index in indices if index in projected]
    if not points:
        return None
    xs = [point[0] for point in points]
    ys = [point[1] for point in points]
    return (
        min(xs) - padding,
        min(ys) - padding,
        max(xs) + padding,
        max(ys) + padding,
    )


def _initial_backbone_color(
    project: TiamatProject,
    selected_indices: set[int],
) -> str:
    strand_colors = {strand.strand_id: strand.color for strand in project.strands}
    for index in sorted(selected_indices):
        if index not in project.bases:
            continue
        return _backbone_color_for_base(project, index, strand_colors)
    return "#cbd5e1"


def _strand_end_color(
    project: TiamatProject,
    strand,
    strand_colors: dict[int, str],
) -> str:
    if len(strand.base_indices) < 2:
        return strand_colors.get(strand.strand_id, "#cbd5e1")
    previous_index = strand.base_indices[-2]
    return _backbone_color_for_base(project, previous_index, strand_colors)


def _display_recent_file_path(path: str, max_length: int = 34) -> str:
    text = str(path)
    if len(text) <= max_length:
        return text
    if max_length <= 15:
        return f"...{text[-(max_length - 3):]}"
    keep = max_length - 3
    head = min(max(10, keep // 2), keep - 12)
    tail = keep - head
    return f"{text[:head]}...{text[-tail:]}"


def _menu_toggle_mark(enabled: bool) -> str:
    return "\u2713" if enabled else ""


def _disabled_toolbar_icon(image: tk.PhotoImage) -> tk.PhotoImage:
    disabled = tk.PhotoImage(width=image.width(), height=image.height())
    for y in range(image.height()):
        for x in range(image.width()):
            if image.transparency_get(x, y):
                disabled.transparency_set(x, y, True)
                continue
            red, green, blue = _photo_rgb(image.get(x, y))
            gray = int(round(0.299 * red + 0.587 * green + 0.114 * blue))
            gray = int(round(gray * 0.45 + 0xE5 * 0.55))
            disabled.put(f"#{gray:02x}{gray:02x}{gray:02x}", (x, y))
    return disabled


def _photo_rgb(color_value) -> tuple[int, int, int]:
    if isinstance(color_value, tuple):
        return tuple(int(channel) for channel in color_value[:3])
    text = str(color_value).strip()
    if text.startswith("#") and len(text) == 7:
        return (int(text[1:3], 16), int(text[3:5], 16), int(text[5:7], 16))
    if text.startswith("#") and len(text) == 13:
        return (int(text[1:5], 16) // 257, int(text[5:9], 16) // 257, int(text[9:13], 16) // 257)
    parts = text.split()
    if len(parts) >= 3:
        return (int(parts[0]), int(parts[1]), int(parts[2]))
    return (0, 0, 0)


def _selection_targets_for_mode(project: TiamatProject, index: int, mode: str) -> set[int]:
    if index not in project.bases:
        return set()
    if mode in {"base", "box"}:
        return {index}
    if mode == "base_pair":
        targets = {index}
        across = project.bases[index].across
        if across is not None and across in project.bases:
            targets.add(across)
        return targets
    if mode == "strand":
        strand_id = project.bases[index].strand_id
        return _strand_indices(project, strand_id)
    if mode == "helix":
        return _helix_indices(project, index)
    return {index}


def _strand_indices(project: TiamatProject, strand_id: int | None) -> set[int]:
    if strand_id is None:
        return set()
    for strand in project.strands:
        if strand.strand_id == strand_id:
            return set(strand.base_indices)
    return set()


def _helix_indices(project: TiamatProject, index: int) -> set[int]:
    start_strand_id = project.bases[index].strand_id
    if start_strand_id is None:
        return {index}
    adjacency: dict[int, set[int]] = {}
    for base in project.bases.values():
        if base.across is None or base.across not in project.bases:
            continue
        left = base.strand_id
        right = project.bases[base.across].strand_id
        if left is None or right is None:
            continue
        adjacency.setdefault(left, set()).add(right)
        adjacency.setdefault(right, set()).add(left)

    strand_ids: set[int] = set()
    stack = [start_strand_id]
    while stack:
        strand_id = stack.pop()
        if strand_id in strand_ids:
            continue
        strand_ids.add(strand_id)
        stack.extend(sorted(adjacency.get(strand_id, set()) - strand_ids, reverse=True))

    targets: set[int] = set()
    for strand_id in strand_ids:
        targets.update(_strand_indices(project, strand_id))
    return targets


def _projected_indices_in_box(
    projected: dict[int, tuple[float, float, float]],
    start: tuple[int, int],
    end: tuple[int, int],
) -> set[int]:
    left = min(start[0], end[0])
    right = max(start[0], end[0])
    top = min(start[1], end[1])
    bottom = max(start[1], end[1])
    return {
        index
        for index, (x, y, _depth) in projected.items()
        if left <= x <= right and top <= y <= bottom
    }


def _nearest_projected_index(
    projected: dict[int, tuple[float, float, float]],
    x: float,
    y: float,
    max_distance_sq: float = 225.0,
) -> int | None:
    if not projected:
        return None
    best_index = None
    best_distance = float("inf")
    for index, (px, py, _depth) in projected.items():
        distance = (px - x) ** 2 + (py - y) ** 2
        if distance < best_distance:
            best_distance = distance
            best_index = index
    if best_index is None or best_distance > max_distance_sq:
        return None
    return best_index


def _base_edit_options(molecule: str) -> tuple[tuple[str, str], ...]:
    if str(molecule).upper() == "RNA":
        return (
            ("Adenine", "A"),
            ("Uracil", "U"),
            ("Cytosine", "C"),
            ("Guanine", "G"),
            ("Generic", "N"),
        )
    return (
        ("Adenine", "A"),
        ("Thymine", "T"),
        ("Cytosine", "C"),
        ("Guanine", "G"),
        ("Generic", "N"),
    )


def _allowed_sequence_bases(molecule: str) -> tuple[str, ...]:
    if str(molecule).upper() == "RNA":
        return ("A", "U", "C", "G", "N")
    return ("A", "T", "C", "G", "N")


def _sequence_legend_items(molecule: str) -> tuple[str, ...]:
    if str(molecule).upper() == "RNA":
        return (
            "A - Adenine",
            "U - Uracil",
            "C - Cytosine",
            "G - Guanine",
            "N - Generic",
        )
    return (
        "A - Adenine",
        "T - Thymine",
        "C - Cytosine",
        "G - Guanine",
        "N - Generic",
    )


def _normalize_sequence_input(text: str, allowed_bases: tuple[str, ...]) -> tuple[str, list[str]]:
    allowed = set(allowed_bases)
    normalized: list[str] = []
    invalid: list[str] = []
    for raw_char in text.upper():
        if raw_char.isspace():
            continue
        if raw_char in allowed:
            normalized.append(raw_char)
            continue
        invalid.append(raw_char)
    return "".join(normalized), invalid


def _selected_strand_for_indices(project: TiamatProject, selected_indices: set[int]):
    selected = sorted(index for index in selected_indices if index in project.bases)
    if not selected:
        return None
    strand_id = project.bases[selected[0]].strand_id
    if strand_id is None:
        return None
    for index in selected[1:]:
        if project.bases[index].strand_id != strand_id:
            return None
    for strand in project.strands:
        if strand.strand_id == strand_id:
            return strand
    return None


def _basepair_dash(project: TiamatProject, left: int, right: int) -> tuple[int, int] | None:
    left_base = project.bases[left]
    right_base = project.bases[right]
    if basepair_state(left_base, right_base) == "wobble":
        return (4, 2)
    return None


def _svg_dash_attr(dash: tuple[int, int] | None) -> str:
    if dash is None:
        return ""
    return f' stroke-dasharray="{dash[0]},{dash[1]}"'


def _clear_selection_on_empty_click(mode: str) -> bool:
    return mode in EMPTY_CLICK_CLEAR_MODES


def _has_shortcut_modifier(state: int) -> bool:
    return bool(state & (0x0004 | 0x0008 | 0x0080))


def _has_shift_modifier(state: int) -> bool:
    return bool(state & 0x0001)


def _selected_edge_sets(
    project: TiamatProject,
    selected_indices: set[int],
) -> tuple[list[tuple[int, int]], list[tuple[int, int]]]:
    if not selected_indices:
        return [], []
    selected_along: list[tuple[int, int]] = []
    for strand in project.strands:
        for left, right in zip(strand.base_indices, strand.base_indices[1:]):
            if left in selected_indices:
                selected_along.append((left, right))

    selected_across: list[tuple[int, int]] = []
    seen_across: set[tuple[int, int]] = set()
    for index in sorted(selected_indices):
        across = project.bases[index].across
        if across is None:
            continue
        edge = tuple(sorted((index, across)))
        if edge in seen_across:
            continue
        seen_across.add(edge)
        selected_across.append(edge)
    return selected_along, selected_across


def _widget_is_descendant(widget: tk.Widget, ancestor: tk.Widget) -> bool:
    current = widget
    while True:
        if current is ancestor:
            return True
        parent_name = current.winfo_parent()
        if not parent_name:
            return False
        parent = current.nametowidget(parent_name)
        if parent is current:
            return False
        current = parent
