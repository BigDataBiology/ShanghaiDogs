from jug import TaskGenerator
from jug.utils import cached_glob


@TaskGenerator
def compress_emapper_file(emapper):
    import subprocess
    from os import path
    import os
    # --threads=4 both speeds up the compression and ensures that the resulting
    # file is chunked (enable parallel decompression)
    subprocess.check_call(['xz', '--threads=4', emapper])

    os.link(emapper + '.xz', '../data/ShanghaiDogsTables/MAGAnnotations/' + path.basename(emapper) + '.xz')
    return emapper + '.xz'


emappers = cached_glob('../intermediate-outputs/eggNOG-annot/*.emapper.annotations')
for emapper in emappers:
    compress_emapper_file(emapper)
