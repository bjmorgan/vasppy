"""Class for storing additional VASP calculation metadata from a YAML file."""

import yaml  # type: ignore


class VASPMeta:
    """Stores additional metadata for a VASP calculation directory."""

    def __init__(
        self,
        title: str,
        description: str,
        status: str,
        notes: str | None = None,
        type: str | None = None,
        track: list[str] | None = None,
    ) -> None:
        """Initialise a VASPMeta object.

        Args:
            title: The title string for this calculation.
            description: Long description of the calculation.
            status: Current status. Must be one of ``'to-run'``,
                ``'incomplete'``, ``'finished'``, or ``'dropped'``.
            notes: Any additional notes. Defaults to None.
            type: Optional calculation type descriptor. Must be one of
                ``'single-point'`` or ``'neb'`` if provided.
                Defaults to None.
            track: Optional list of filenames whose md5 checksums are
                computed when summarising the calculation output.
                Defaults to None.

        Raises:
            ValueError: If *status* is not one of the expected values.
            ValueError: If *type* is provided but not one of the expected values.
        """
        self.title = title
        self.description = description
        self.notes = notes
        expected_status = ["to-run", "incomplete", "finished", "dropped"]
        if status not in expected_status:
            raise ValueError(
                f'Unexpected calculation status: "{status}"'
                f" for calculation {title}"
            )
        self.status = status
        expected_types = ["single-point", "neb"]
        if type:
            if type not in expected_types:
                raise ValueError(
                    f'Unexpected calculation type: "{type}"'
                    f" for calculation {title}"
                )
            self.type: str | None = type
        else:
            self.type = None
        self.track = track

    @classmethod
    def from_file(cls, filename: str) -> "VASPMeta":
        """Create a VASPMeta object by reading a ``vaspmeta.yaml`` file.

        Args:
            filename: Path of the YAML file to read.

        Returns:
            The populated VASPMeta object.
        """
        with open(filename, "r") as stream:
            data = yaml.load(stream, Loader=yaml.SafeLoader)
            notes = data.get("notes")
            v_type = data.get("type")
            track = data.get("track")
            xargs: dict = {}
            if track:
                if isinstance(track, str):
                    track = [track]
                xargs["track"] = track
            vaspmeta = VASPMeta(
                data["title"],
                data["description"],
                data["status"],
                notes=notes,
                type=v_type,
                **xargs,
            )
        return vaspmeta
