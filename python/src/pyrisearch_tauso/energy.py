"""Grouped energy statistics over score-filtered RIsearch hits."""

from typing import NamedTuple

import pyarrow as pa
import pyarrow.compute as pc

from . import Reduction


class EnergyStats(NamedTuple):
    """Boltzmann weight sum and minimum reported energy, in kcal/mol."""

    sum_exp: float
    min_energy: float


def energy_stats(cutoffs, *, group_by=("query", "target"), rt):
    """Reduce hits to {cutoff: {key: EnergyStats}} using score > cutoff.

    rt is the caller's RT in kcal/mol. group_by selects query, target, or both;
    a single column gives scalar keys, multiple columns give tuple keys.
    Sums are computed per batch then summed across batches. Like other floating
    sums, their last bits can depend on batch boundaries and Arrow threading.
    No hits returns {}. A cutoff without qualifying hits has an empty mapping
    when other hits were read. Filtering several cutoffs from one search requires
    min_score <= min(cutoffs); use neighborhood=0 for independent endpoint hits.
    """
    cutoffs = sorted({int(c) for c in cutoffs})
    keys = tuple(group_by)
    if not keys or len(set(keys)) != len(keys) or any(k not in ("query", "target") for k in keys):
        raise ValueError("group_by must contain query and/or target without duplicates")
    if not 0 < rt < float("inf"):
        raise ValueError("rt must be finite and positive")

    schema = {
        "cutoff": pa.int64(),
        **{k: pa.string() for k in keys},
        "sum_exp": pa.float64(),
        "min_energy": pa.float64(),
    }

    def combine(batch):
        base = pa.table(
            {
                **{k: batch.column(k) for k in keys},
                "score": batch.column("score"),
                "_exp": pc.exp(pc.divide(batch.column("energy"), -rt)),
                "energy": batch.column("energy"),
            }
        )
        parts = []
        for cutoff in cutoffs:
            sub = base.filter(pc.greater(base.column("score"), cutoff))
            if not sub.num_rows:
                continue
            agg = sub.group_by(list(keys)).aggregate([("_exp", "sum"), ("energy", "min")])
            parts.append(
                pa.table(
                    {
                        "cutoff": pa.array([cutoff] * agg.num_rows, pa.int64()),
                        **{k: agg.column(k) for k in keys},
                        "sum_exp": agg.column("_exp_sum"),
                        "min_energy": agg.column("energy_min"),
                    }
                )
            )
        return (
            pa.concat_tables(parts)
            if parts
            else pa.table({name: pa.array([], dtype) for name, dtype in schema.items()})
        )

    def finalize(table):
        merged = table.group_by(["cutoff", *keys]).aggregate([("sum_exp", "sum"), ("min_energy", "min")])
        result = {c: {} for c in cutoffs}
        columns = [merged.column(k).to_pylist() for k in ("cutoff", *keys, "sum_exp_sum", "min_energy_min")]
        for cutoff, *values in zip(*columns):
            key = values[0] if len(keys) == 1 else tuple(values[:-2])
            result[cutoff][key] = EnergyStats(float(values[-2]), float(values[-1]))
        return result

    return Reduction(columns=(*keys, "score", "energy"), combine=combine, finalize=finalize, empty={})
