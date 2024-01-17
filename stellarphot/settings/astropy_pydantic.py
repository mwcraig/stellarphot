from typing import Annotated

from astropy.units import IrreducibleUnit, Quantity, Unit
from pydantic import PlainSerializer, PlainValidator, WithJsonSchema
from typing_extensions import TypeAliasType

__all__ = ["UnitType", "QuantityType", "PixelScaleType"]

# Approach to validation of units was inspired by the GammaPy project
# which did it before we did:
# https://docs.gammapy.org/dev/_modules/gammapy/analysis/config.html


def ser(v):
    return str(v)


def schema_thing(value, handler, **kwd):
    # print(f"type(value)={type(value)}")
    print(f"{value=}")
    print(f"{handler=}")
    print(f"{kwd=}")
    # print(f"{value=}")
    # print(dir(handler))
    # reals = [h for h in dir(handler) if not h.startswith("_")]
    # for real in reals:
    #     print(f"{real=}")
    #     print([h for h in dir(getattr(handler, real)) if not h.startswith("_")])
    return {"type": "string"}


UnitType = TypeAliasType(
    "UnitType",
    Annotated[
        Unit | IrreducibleUnit,
        # Must use PlainValidator below to completely replace pydantic validation,
        # because pydantic has no idea how to validate a Unit, and the other validators
        # invoke the pydantic validation at some point.
        PlainValidator(lambda v: Unit(v)),
        # GetPydanticSchema(get_pydantic_json_schema=schema_thing),
        # # You need some kind of schema for pydantic to process this annotation
        WithJsonSchema({"type": "string", "examples": ["m", "m / s"]}),
        # Serialize to a string when dumping to JSON
        PlainSerializer(ser, when_used="json"),
    ],
)
# class UnitType(Unit):
#     # Validator for Unit type
#     @classmethod
#     def __get_validators__(cls):
#         yield cls.validate

#     @classmethod
#     def validate(cls, v):
#         return Unit(v)

#     @classmethod
#     def __modify_schema__(cls, field_schema, field):
#         # Set default values for the schema in case the field doesn't provide them
#         name = "Unit"
#         description = "An astropy unit"

#         name = field.name or name
#         description = field.field_info.description or description
#         examples = field.field_info.extra.get("examples", [])

#         field_schema.update(
#             {
#                 "title": name,
#                 "description": description,
#                 "examples": examples,
#                 "type": "string",
#             }
#         )


class QuantityType(Quantity):
    # Validator for Quantity type
    @classmethod
    def __get_validators__(cls):
        yield cls.validate

    @classmethod
    def validate(cls, v):
        try:
            v = Quantity(v)
        except TypeError as err:
            raise ValueError(f"Invalid value for Quantity: {v}") from err
        else:
            if not v.unit.bases:
                raise ValueError("Must provided a unit")
        return v

    # @classmethod
    # def __modify_schema__(cls, field_schema, field):
    #     # Set default values for the schema in case the field doesn't provide them
    #     name = "Quantity"
    #     description = "An astropy Quantity with units"

    #     name = field.name or name
    #     description = field.field_info.description or description
    #     examples = field.field_info.extra.get("examples", [])

    #     field_schema.update(
    #         {
    #             "title": name,
    #             "description": description,
    #             "examples": examples,
    #             "type": "string",
    #         }
    #     )


class PixelScaleType(Quantity):
    # Validator for pixel scale type
    @classmethod
    def __get_validators__(cls):
        yield cls.validate

    @classmethod
    def validate(cls, v):
        try:
            v = Quantity(v)
        except TypeError as err:
            raise ValueError(f"Invalid value for Quantity: {v}") from err
        if (
            len(v.unit.bases) != 2
            or v.unit.bases[0].physical_type != "angle"
            or v.unit.bases[1].name != "pix"
        ):
            raise ValueError(f"Invalid unit for pixel scale: {v.unit!r}")
        return v

    # @classmethod
    # def __modify_schema__(cls, field_schema):
    #     field_schema.update(
    #         {
    #             "title": "PixelScale",
    #             "description": "An astropy Quantity with units of angle per pixel",
    #             "examples": ["0.563 arcsec / pix"],
    #             "type": "string",
    #         }
    #     )
