import imageio
from glob import iglob


def animate(output, format, pause):
    """
    Create animation from pack of pictures.
    For our knowledge only png format amongst widely spread is supported
    :param output: str - path to directory with images
    :param format: str - format of using pictures
    :param pause: float = intervals in seconds for frames
    :return:
    """
    appropriate_formats = ['png']
    if format not in appropriate_formats:
        format = 'png'
    with imageio.get_writer(f'{output}/assembly.gif', mode='I', duration=pause) as writer:
        for image in sorted(iglob(f'{output}/*.{format}')):
            frame = imageio.imread(image)
            writer.append_data(frame)

