from chimerax.core.toolshed import BundleAPI


class _MyAPI(BundleAPI):
    @staticmethod
    def start_tool(session, tool_name):
        from .tool import StrudelScore
        # return StrudelScore(session, tool_name)
        # return tool.StrudelScore(session, tool_name)
        from chimerax.core import tools
        return tools.get_singleton(session, StrudelScore, tool_name, create=True)



    @staticmethod
    def register_command(command_name, logger):
        from . import cmd
        if command_name == 'strudel':
            from . import cmd
            cmd.register_strudel(logger)



bundle_api = _MyAPI()

