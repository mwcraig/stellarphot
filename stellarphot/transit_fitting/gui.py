import json
import re
from pathlib import Path

import ipywidgets as ipw
import numpy as np
from astropy.utils.data import get_pkg_data_filename
from traitlets import Bool, observe

from stellarphot.transit_fitting.io import get_tic_info

__all__ = [
    "MyValid",
    "make_checker",
    "validate_exposure_time",
    "populate_TIC_boxes",
    "populate_TOI_boxes",
    "exotic_settings_widget",
    "set_values_from_json_file",
    "get_values_from_widget",
    "generate_json_file_name",
]


template_types = ["known", "candidate"]
template_json = {}
to_fill = {}

for template in template_types:
    template_name = get_pkg_data_filename(
        "data/tic-template-for-exotic-" f"{template}.json"
    )
    with open(template_name) as f:
        template_json[template] = json.load(f)

    template_name = get_pkg_data_filename(f"data/exotic-to-mod-{template}.json")
    with open(Path(template_name)) as f:
        to_fill[template] = json.load(f)

exotic_arguments = dict(
    known=["--nasaexoarch", "--pre"], candidate=["--override", "--pre"]
)

# Nested keys are flattened by joining them with this character
join_char = "😬"


class MyValid(ipw.Button):
    """
    A class containing a more compact indicator of valid entries based on ipywidgets
    button value.

    Parameters
    ----------

    **kwd :
        All keyword arguments used to initialize an instance of this class are
        passed to the ipywidgets.Button constructor.

    Attributes
    ----------

    value : bool
        Current value of the indicator. Initialized to False.

    """

    value = Bool(False, help="Bool value").tag(sync=True)

    def __init__(self, **kwd):
        super().__init__(**kwd)
        self.layout.width = "40px"
        self._set_properties(None)

    @observe("value")
    def _set_properties(self, change):  # noqa: ARG002
        """
        Widget callbacks need to accept a single argument, even if it is not used.
        """
        if self.value:
            self.style.button_color = "green"
            self.icon = "check"
        else:
            self.style.button_color = "red"
            self.icon = "times"


def make_checker(indicator_widget, value_widget):
    """
    Build an observer that checks TIC number and, if it is a valid
    TIC number, looks up information about the TIC object from
    MAST copy of TIC as priors for EXOTIC. It also sets the visible checkbox
    to the appropriate state.

    Parameters
    ----------

    indicator_widget : `~stellarphot.transit_fitting.gui.MyValid` widget
        The widget that indicates to the user whether or not the value is
        reasonable.

    value_widget: ipywidget
        This widget should be generated by exotic_settings_widget.

    Returns
    -------

    function
        Function with the correct signature for use as an observer on an
        ipywidget.
    """

    def check_name(change):
        # Valid TIC number is 9 digits
        ticced = re.compile(r"TIC \d{9,10}$")
        owner = change["owner"]
        is_tic = ticced.match(change["new"])
        if is_tic:
            if indicator_widget is not None:
                indicator_widget.value = True
            owner.disabled = True
            tic_info = get_tic_info(change["new"][-9:])
            if not tic_info:
                indicator_widget.value = False
                indicator_widget.tooltip = "Not a valid TIC number"
            else:
                populate_TIC_boxes(tic_info, value_widget)
            owner.disabled = False
        else:
            owner.disabled = False
            if indicator_widget is not None:
                indicator_widget.value = False
                indicator_widget.tooltip = "TIC numbers have 9 digits"

    return check_name


def validate_exposure_time(
    indicator_widget, value_widget  # noqa: ARG001 (value_widget needed for callback)
):
    """Validates the exposure time input.

    Parameters
    ----------
    indicator_widget : `~stellarphot.transit_fitting.gui.MyValid` widget
        The widget that indicates to the user whether or not the value is
        reasonable.

    value_widget: ipywidget
        This widget should be generated by exotic_settings_widget.

    Returns
    -------
    function
        Function that will set the correct boolean value on the
        indicator_widget to indicate if the value of the exposure time
        is valid.  This can be used as an observer for an ipywidget
    """

    def check_exposure(change):
        # Valid Exposure time is greater than zero
        if change["new"] > 0:
            if indicator_widget is not None:
                indicator_widget.value = True
        else:
            if indicator_widget is not None:
                indicator_widget.value = False

    return check_exposure


def populate_TIC_boxes(tic_info, value_widget):
    """
    Set the appropriate widget values given information pulled from
    MAST TIC.

    Parameters
    ----------

    tic_info : dict
        Dictionary of information about the TIC object.

    value_widget: ipywidget
        This widget should contain a 'candidate' key.

    Returns
    -------

    None
        Sets values of planetary parameters of `candidate` in ``value_widget`` in place.
    """
    # Match EXOTIC json keys to columns in the result returned from
    # astroquery
    exotic_tic = {
        "Star Effective Temperature (K)": "Teff",
        "Star Effective Temperature (+) Uncertainty": "epos_Teff",
        "Star Effective Temperature (-) Uncertainty": "eneg_Teff",
        "Star Surface Gravity (log(g))": "logg",
        "Star Surface Gravity (+) Uncertainty": "epos_logg",
        "Star Surface Gravity (-) Uncertainty": "eneg_logg",
        "Host Star Name": "UCAC",
        "Star Metallicity ([FE/H])": "MH",
        "Star Metallicity (+) Uncertainty": "e_MH",
        "Star Metallicity (-) Uncertainty": "e_MH",
    }
    for k, v in exotic_tic.items():
        exotic_key = join_char.join(["planetary_parameters", k])
        if k == "Host Star Name":
            value_widget["candidate"][exotic_key].value = f"UCAC4 {tic_info[v][0]}"
        elif not np.isnan(tic_info[v][0]):
            value_widget["candidate"][exotic_key].value = tic_info[v][0]


def populate_TOI_boxes(toi, exotic_widget):
    """
    Set the appropriate widget values given information pulled from
    TESS Target of Interest (TOI) list.

    Parameters
    ----------

    toi : `~stellarphot.io.tess.TOI`
        Information about the TESS Target of Interest (TOI) object.

    exotic_widget:  ipywidget
        This widget should be generated by exotic_settings_widget.

    Returns
    -------

    None
        Sets values of planetary parameters of `candidate` in
        ``exotic_widget`` in place.
    """
    # Match EXOTIC json keys to columns in the result returned from
    # astroquery
    exotic_toi = {
        "Planet Name": "tic_id",
        "Target Star RA": "coord",
        "Target Star Dec": "coordP",
        "Orbital Period (days)": "period",
        "Orbital Period Uncertainty": "period_error",
        "Published Mid-Transit Time (BJD-UTC)": "epoch",
        "Mid-Transit Time Uncertainty": "epoch_error",
        # Could maybe get these from TOI information, but not straightforward
        # "Ratio of Planet to Stellar Radius (Rp/Rs)": 0.0,
        # "Ratio of Planet to Stellar Radius (Rp/Rs) Uncertainty": 0.0,
        # "Ratio of Distance to Stellar Radius (a/Rs)": 0.0,
        # "Ratio of Distance to Stellar Radius (a/Rs) Uncertainty": 0.0,
    }
    for k, v in exotic_toi.items():
        exotic_key = join_char.join(["planetary_parameters", k])
        if k == "Planet Name":
            exotic_widget["candidate"][exotic_key].value = f"TIC {toi.tic_id}"
        elif k == "Target Star RA":
            exotic_widget["candidate"][exotic_key].value = toi.coord.ra.to_string(
                unit="hour", decimal=False, sep=":"
            )
        elif k == "Target Star Dec":
            exotic_widget["candidate"][exotic_key].value = toi.coord.dec.to_string(
                unit="degree", decimal=False, sep=":"
            )
        else:
            exotic_widget["candidate"][exotic_key].value = getattr(toi, v).value


# This sets up the specific widgets whose values get validated.
# That includes the widget that contains the TIC number for candidates.
validators = dict(known={}, candidate={})
validators["candidate"]["Planet Name"] = make_checker
for k in validators:
    validators[k]["Exposure Time (s)"] = validate_exposure_time


def exotic_settings_widget():
    """
    Generate a widget to enter (and store) settings for exotic.

    Parameters
    ----------

    None

    Returns
    -------

    ipywidget
        Widget with settings for EXOTIC.
    """

    # We rely on some global variables:
    global to_fill, template_types

    # This dictionary will contain all of the widgets
    widget_list = {}

    # Each widget has the same layout for its description and its value/input
    layout_description = ipw.Layout(width="70%")
    layout_input = ipw.Layout(width="30%")

    # Maintain a separate dict of just the value widgets
    value_widget = {}

    # For exotic there are two templates, one for known exoplanets and one for
    # candidate exoplanets
    for template in template_types:
        # Keep separate track of sidgets that make up each template
        value_widget[template] = {}
        widget_list[template] = []
        # Autogenerate the input fields and labels
        for k in to_fill[template]:
            for k2, v in to_fill[template][k].items():
                if isinstance(v, str):
                    input_widget = ipw.Text(value=v)
                elif isinstance(v, float):
                    if v >= 0:
                        # In principle we could limit this to positive values
                        input_widget = ipw.FloatText(value=v)
                    else:
                        input_widget = ipw.FloatText(value=v)
                input_widget.layout = layout_input

                # Each horizontal box below is one "cell" and the input grid.
                # The HTML widget has the label, and input_widget is the text or
                # float widget.
                hb = ipw.HBox(
                    [ipw.HTML(value=k2, layout=layout_description), input_widget]
                )

                # Widgets with validators need some observers added
                try:
                    # validator is a function that generates a function with the right
                    # with the right signature to handle an ipywidgets event
                    validator = validators[template][k2]
                except KeyError:
                    # Most of the widgets...
                    pass
                else:
                    kids = list(hb.children)
                    # Add a visible indication of whether the value is reasonable
                    indicator = MyValid(value=False)
                    kids.append(indicator)

                    # Add an observer to the widget to watch for changes
                    input_widget.observe(
                        validator(indicator, value_widget), names="value"
                    )

                    hb.children = kids

                val_key = join_char.join([k, k2])
                value_widget[template][val_key] = input_widget
                widget_list[template].append(hb)

    hb2 = {}
    for template in template_types:
        hb2[template] = ipw.HBox(
            [
                ipw.VBox(widget_list[template][:16], layout=ipw.Layout(padding="10px")),
                ipw.VBox(widget_list[template][16:]),
            ]
        )

    select_planet_type = ipw.ToggleButtons(
        description="Known or candidate exoplanet?",
        options=template_types,
        style={"description_width": "initial"},
    )

    lookup_link_text = dict(
        known="https://exoplanetarchive.ipac.caltech.edu/",
        candidate="https://exofop.ipac.caltech.edu/tess/",
    )

    lookup_link_html = {}

    for k, v in lookup_link_text.items():
        lookup_link_html[k] = ipw.HTML(
            f"<h3>For some information about this "
            f'object: <a href="{v}" target="_blank">{v}</a></h3>'
        )

    input_container = ipw.VBox()

    whole_thing = ipw.VBox(children=[select_planet_type, input_container])
    whole_thing.planet_type = select_planet_type
    whole_thing.value_widget = value_widget
    pre_reduced_file = join_char.join(["optional_info", "Pre-reduced File:"])
    whole_thing.data_file_widget = {
        "candidate": value_widget["candidate"][pre_reduced_file],
        "known": value_widget["known"][pre_reduced_file],
    }

    def observe_select(change):  # noqa: ARG001
        """
        Widget callbacks need to accept a single argument, even if it is not used.
        """
        input_container.children = [
            lookup_link_html[select_planet_type.value],
            hb2[select_planet_type.value],
        ]

    select_planet_type.observe(observe_select, names="value")
    observe_select(select_planet_type.value)

    return whole_thing


def set_values_from_json_file(widget, json_file):
    """
    Set values in EXOTIC widget from a JSON file of those settings.

    Parameters
    ----------

    widget : ipywidget
        This widget should be generated by exotic_settings_widget.

    json_file : str
        File with settings for the widget.

    Returns
    -------
    None
        Sets values of parameters of widget in place.
    """
    with open(json_file) as f:
        input_values = json.load(f)

    planet_type = widget.planet_type.value
    for k, a_widget in widget.value_widget[planet_type].items():
        k1, k2 = k.split(join_char)
        a_widget.value = input_values[k1][k2]


def get_values_from_widget(exotic_widget, key=None):
    """
    Extract EXOTIC settings from widget.

    Parameters
    ----------

    exotic_widget : ipywidget
        This widget should be generated by exotic_settings_widget.

    key : str, either "known" or "candidate", optional
        Indicates which case to use for EXOTIC. If ``None``, use the ``planet_type``
        attribute of the ``exotic_widget``.

    Returns
    -------

    dict
        Value of modified template_json[key], where template_json is a global variable.

    """
    if not key:
        key = exotic_widget.planet_type.value

    for k, widget in exotic_widget.value_widget[key].items():
        k1, k2 = k.split(join_char)
        template_json[key][k1][k2] = widget.value

    return template_json[key]


def generate_json_file_name(exotic_widget, key=None):
    """
    Generate settings file name from user input.

    Parameters
    ----------

    exotic_widget : ipywidget
        This widget should be generated by exotic_settings_widget.

    key : str, either "known" or "candidate", optional
        Indicates which case to use for EXOTIC. If ``None``, use the ``planet_type``
        attribute of the ``exotic_widget``.

    Returns
    -------
    str
        Name of file to save settings to.
    """
    if not key:
        key = exotic_widget.planet_type.value

    get_values_from_widget(exotic_widget, key=key)
    user_info = "user_info"
    planet = "planetary_parameters"
    filter_key = "Filter Name (aavso.org/filters)"
    date = template_json[key][user_info]["Observation date"]
    planet_name = template_json[key][planet]["Planet Name"]
    filter_name = template_json[key][user_info][filter_key]
    name = f"{planet_name}-{date}-{filter_name}"
    return name.replace(" ", "_") + ".json"
