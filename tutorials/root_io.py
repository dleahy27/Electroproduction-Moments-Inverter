"""Merge scan point trees into one file using RDataFrame snapshots.

The original point files remain available for fit inputs and resuming runs.
A scan_case column keeps every merged row tied to its manifest entry.
"""

from array import array
from pathlib import Path
from tempfile import TemporaryDirectory
import ROOT
from tutorials.analysis_utils import truth_parameters
from tutorials.paths import PROJECT_DIR


def merge_scan(points, output, trees, mass_truth_key=None):
    """Write the manifest's trees and scan coordinates to one atomic ROOT file.

    trees maps each tree name to its source-path key in a manifest row. ROOT
    handles tree concatenation; no manual branch-address copying is needed.
    mass_truth_key optionally copies scalar BW/background truth coordinates
    into a massModelTruth tree so plotting never reopens individual truth files.
    """
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary_output = output.with_name(output.stem + ".merging.root")
    with TemporaryDirectory(prefix=".merge_", dir=output.parent) as temporary:
        merger = ROOT.TFileMerger(False, False)
        merger.OutputFile(str(temporary_output), "RECREATE")
        for case, point in enumerate(points):
            part = Path(temporary) / f"case_{case:04d}.root"
            coordinates = {"scan_case": int(point.get("scan_case", case))}
            for field in (
                "mass_index",
                "asymmetry_index",
                "mass_GeV",
                "k_asymmetry",
                "repeat",
                "suppression",
            ):
                if field in point:
                    coordinates["scan_" + field] = float(point[field])
            for tree_name, key in trees.items():
                source = PROJECT_DIR / point[key]
                frame = ROOT.RDataFrame(tree_name, str(source))
                for name, value in coordinates.items():
                    frame = frame.Define(name, repr(value))
                options = ROOT.RDF.RSnapshotOptions()
                options.fMode = "UPDATE" if part.exists() else "RECREATE"
                frame.Snapshot(tree_name, str(part), "", options)
            if mass_truth_key:
                parameters = truth_parameters(PROJECT_DIR / point[mass_truth_key])
                destination = ROOT.TFile.Open(str(part), "UPDATE")
                truth = ROOT.TTree("massModelTruth", "Injected mass-model coordinates")
                storage = {
                    name: array("d", [value])
                    for name, value in {**parameters, **coordinates}.items()
                }
                for name, value in storage.items():
                    truth.Branch(name, value, name + "/D")
                truth.Fill()
                truth.Write()
                destination.Close()
            if not merger.AddFile(str(part)):
                raise OSError(f"Cannot merge scan case {case}")
        if not merger.Merge():
            raise OSError(f"Cannot merge scan into {output}")
    temporary_output.replace(output)
