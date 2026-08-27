#  This source file is part of the Avogadro project
#  This source code is released under the 3-Clause BSD License, (see "LICENSE").
#  https://github.com/ghutchis/avogadro-stk/

"""Dispatch a feature identifier to its handler.

Every path out of here is JSON. A traceback on stdout would leave Avogadro
with nothing to show but a parse error, so anything that escapes a feature is
turned into an `error` the user can act on; `--debug` appends the traceback.
"""

import json
import sys
import traceback


def dispatch(feature, avo_input, options):
    match feature:
        case "analyze":
            from .features.analyze import run as handler
            return handler(avo_input)
        case "polymer":
            from .features.chains import run as handler
            return handler(avo_input, options, cyclic=False)
        case "macrocycle":
            from .features.chains import run as handler
            return handler(avo_input, options, cyclic=True)
        case "cage":
            from .features.cage import run as handler
            return handler(avo_input, options)
        case "cof":
            from .features.cof import run as handler
            return handler(avo_input, options)
    return None


def run(args):
    avo_input = json.load(sys.stdin)
    # Avogadro nests the dialog answers, but older builds have been known to
    # merge them into the top level; accept either.
    options = avo_input.get("options")
    if not isinstance(options, dict):
        options = avo_input

    try:
        output = dispatch(args.feature, avo_input, options)
    except Exception as exc:  # noqa: BLE001 - nothing may escape as a traceback
        message = str(exc) or exc.__class__.__name__
        if args.debug:
            message += "\n\n" + traceback.format_exc()
        output = {"error": message}

    if output is not None:
        print(json.dumps(output))
