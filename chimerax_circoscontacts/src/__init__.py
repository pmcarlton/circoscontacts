from chimerax.core.toolshed import BundleAPI


class _CircosContactsBundleAPI(BundleAPI):
    @staticmethod
    def register_command(command_name, logger):
        from .cmd import register_command

        register_command(command_name, logger)


bundle_api = _CircosContactsBundleAPI()
