"""TEST ONLY: file-backed collective API used when no MPI runtime is installed.

It exercises the complete writer in independent Python processes, including
barriers, halo object factories, metadata gather, and final file merge.  It is
NOT an MPI implementation/transport test and is never imported by production.
"""
import os
import pickle
from pathlib import Path
import runpy
import sys
import time
import types

writer, input_file, rank, size, shared = sys.argv[1:]
rank, size = int(rank), int(size)
shared = Path(shared)
shared.mkdir(parents=True, exist_ok=True)


class TestCommunicator:
    sequence = 0

    def Get_rank(self):
        return rank

    def Get_size(self):
        return size

    def _exchange(self, value):
        stage = self.sequence
        self.sequence += 1
        path = shared / ('%03d_%03d.pkl' % (stage, rank))
        temporary = path.with_suffix('.tmp')
        temporary.write_bytes(pickle.dumps(value))
        os.replace(temporary, path)
        expected = [shared / ('%03d_%03d.pkl' % (stage, r)) for r in range(size)]
        deadline = time.monotonic() + 90
        while not all(p.exists() for p in expected):
            if (shared / 'abort').exists() or time.monotonic() > deadline:
                raise RuntimeError('Test communicator timeout/abort at stage %d' % stage)
            time.sleep(0.01)
        return [pickle.loads(p.read_bytes()) for p in expected]

    def Barrier(self):
        self._exchange(None)

    def allgather(self, value):
        return self._exchange(value)

    def gather(self, value, root=0):
        values = self._exchange(value)
        return values if rank == root else None

    def allreduce(self, value, op=None):
        values = self._exchange(value)
        return max(values) if op == 'MAX' else sum(values)

    def bcast(self, value, root=0):
        return self._exchange(value)[root]

    def Abort(self, code=1):
        (shared / 'abort').write_text(str(code))
        raise SystemExit(code)


module = types.ModuleType('mpi4py')
module.MPI = types.SimpleNamespace(COMM_WORLD=TestCommunicator(), SUM='SUM', MAX='MAX')
sys.modules['mpi4py'] = module
sys.path[:0] = [str(Path(writer).parent), os.getcwd()]
sys.argv = [writer, input_file]
runpy.run_path(writer, run_name='__main__')
