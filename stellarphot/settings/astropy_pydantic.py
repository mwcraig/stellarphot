from typing import Annotated

from astropy.units import IrreducibleUnit, Quantity, Unit
from pydantic import PlainSerializer, PlainValidator, WithJsonSchema
from typing_extensions import TypeAliasType

__all__ = ["UnitType", "QuantityType", "PixelScaleType"]

# Approach to validation of units was inspired by the GammaPy project
# which did it before we did:
# https://docs.gammapy.org/dev/_modules/gammapy/analysis/config.html


def _quantity_unit_serialization(v):
    """
    Serialize a Quantity or Unit to a string.
    """
    return str(v)


# MAKE A CLOSURE HERE TO HANDLE QUANTITY OR UNIT
def _quantity_unit_type_validator(value, info):
    """
    Validate a Quantity or Unit.

    Parameters
    ----------
    value : str
        The value to validate.
    info : `pydantic.FieldInfo`
        Information about the field being validated.
    """
    print(f"{info=}")
    print(f"{value=}")

    # It is import to TRY QUANTITY FIRST because an expression like "3 meter" will
    # be parsed as a CompositeUnit with value "3 m". f you try to parse "meter" as a
    # Quantity a Type Error will be raised.
    try:
        return_val = Quantity(value)
    except TypeError:
        try:
            # Oh foo. An empty string is a valid Unit. Bork that!
            if value == "":
                raise ValueError("Empty string is not a valid unit")
            return_val = Unit(value)
        except TypeError as erru:
            raise ValueError(
                f"Unable to interpret {value} as an astropy Qunatity or Unit"
            ) from erru

    if isinstance(return_val, Quantity) and not return_val.unit.bases:
        raise ValueError(
            "A unit is required. Use the unit 'dimensionless_unscaled' if "
            "no units are desired."
        )

    return return_val


def _pixel_scale_unit_checker(value):
    """
    Check that the unit is an angle per pixel.
    """
    if (
        len(value.unit.bases) != 2
        or value.unit.bases[0].physical_type != "angle"
        or value.unit.bases[1].name != "pix"
    ):
        raise ValueError(f"Invalid unit for pixel scale: {value.unit!r}")
    return value


UnitType = TypeAliasType(
    "UnitType",
    Annotated[
        # Some units, like meter, are reducible, and can be expressed in terms of
        # other units. Some units cannot be reduced. We need to accommodate both.
        Unit | IrreducibleUnit,
        # Must use PlainValidator below to completely replace pydantic validation,
        # because pydantic has no idea how to validate a Unit, and the other validators
        # invoke the pydantic validation at some point.
        PlainValidator(_quantity_unit_type_validator),
        # Serialize to a string when dumping to JSON
        PlainSerializer(_quantity_unit_serialization, when_used="json"),
        # You need some kind of schema for pydantic to process this annotation -- in
        # particular, it needs to know the type of the field as it is encoded in JSON.
        # TODO: Look at what WithJsonSchema does to see if we can get field_name
        # TODO: Look at whether a BeforeValidator might work to get a title JSON schema
        #       NO DOES NOT WORK
        # TODO: Look at documentation about customizing the schema
        WithJsonSchema({"type": "string"}),
    ],
)


# See the UnitType above for a more detailed explanation of the choices below.
QuantityType = TypeAliasType(
    "QuantityType",
    Annotated[
        Quantity,
        PlainValidator(_quantity_unit_type_validator),
        WithJsonSchema({"type": "string"}),
        PlainSerializer(_quantity_unit_serialization, when_used="json"),
    ],
)


PixelScaleType = TypeAliasType(
    "PixelScaleType",
    Annotated[
        QuantityType,
        PlainValidator(_quantity_unit_type_validator),
        WithJsonSchema({"type": "string"}),
        PlainSerializer(_quantity_unit_serialization, when_used="json"),
    ],
)
