"""AutoDock Vina / smina docking fallback — separate from the mCSM-lig client."""


class VinaError(RuntimeError):
    """Vina/smina could not produce a WT-vs-mutant docking delta for this candidate."""


class VinaDockClient:
    """Receptor/ligand prep and WT-vs-mutant docking.

    Structure preparation (protonation, grid box, mutant modeling) is a
    substantial pipeline of its own. This client exists so mCSM-lig fallbacks
    have a distinct, testable call site; it raises VinaError until that work
    is implemented as a dedicated step.
    """

    def score(
        self,
        *,
        pdb_text: str,
        smiles: str,
        mutation: str,
        position: int,
    ):
        raise VinaError(
            "AutoDock Vina/smina WT-vs-mutant docking is not implemented yet "
            "(receptor prep, grid box, and mutant modeling need their own prompt)"
        )
