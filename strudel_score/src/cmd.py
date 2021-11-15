from chimerax.core.commands import CmdDesc, register
from chimerax.core.commands import StringArg


def get_singleton(session, create=True, display=True):
    if not session.ui.is_gui:
        return None
    from chimerax.core import tools
    from .tool import StrudelScore
    return tools.get_singleton(session, StrudelScore, 'Strudel Score',
                               create=create, display=display)


def strudel_start(session, create=True, display=True):
    ''' Start the Strudel GUI '''
    get_singleton(session, create=create, display=display)
    return session.strudel


def strudel_open(session, path):
    strudel_start(session)
    session.strudel.open_project(path)


def set_library(session, path):
    session_elem = dir(session)
    if 'strudel' not in session_elem:
        get_singleton(session, display=False)
    try:
        session.strudel.set_library(path)
    except AttributeError:
        # This will happen when ran with --nogui option
        from .functions import set_library
        set_library(path)

def register_strudel_start(logger):
    desc = CmdDesc(
        synopsis='Start the Strudel GUI'
    )
    register('strudel start', desc, strudel_start, logger=logger)


def register_strudel_open(logger):
    desc = CmdDesc(synopsis='Open strudel validation results',
                   required=[("path", StringArg)])
    register('strudel open', desc, strudel_open, logger=logger)


def register_set_library(logger):
    desc = CmdDesc(synopsis='Set strudel libraries path',
                   required=[("path", StringArg)])
    register('strudel setLib', desc, set_library, logger=logger)


def register_strudel(logger):
    register_strudel_start(logger)
    register_strudel_open(logger)
    register_set_library(logger)