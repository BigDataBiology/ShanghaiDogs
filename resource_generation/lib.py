from contextlib import contextmanager

@contextmanager
def xz_out(fname):
    import subprocess
    with open(fname, 'wb') as f_raw:
        with subprocess.Popen(['xz', '--threads=12'], stdin=subprocess.PIPE, stdout=f_raw) as p:
            try:
                yield p.stdin
            finally:
                p.stdin.close()
                p.wait()


def pad9(prefix, n):
    n = f'{n:09}'
    n = f'{n[:3]}_{n[3:6]}_{n[6:]}'
    return f'{prefix}.{n}'

