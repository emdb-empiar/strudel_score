from chimerax.core.toolshed import BundleAPI


class _MyAPI(BundleAPI):
    @staticmethod
    def start_tool(session, tool_name):
        from .tool import StrudelScore
        return StrudelScore(session, tool_name)


bundle_api = _MyAPI()
