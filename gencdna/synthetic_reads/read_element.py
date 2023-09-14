"""Class for creating read elements.

Read elements refer to introns, exons or their concatenated combinations.
"""

from typing import Self, TypeAlias


class ReadElement(object):

    def __init__(
        self,
        label: str = '',
        sequence: str = '',
    ) -> None:
        self.label = label
        self.sequence = sequence
        self.length = len(sequence)

    def to_fasta_str(self) -> str:
        return f'>{self.label}\n{self.sequence}\n'

    def __add__(self, other: TypeAlias = Self) -> TypeAlias:
        """Add two ReadElements by concatenating sequences and adding lengths.

        For example:
        a = ReadElement('seq1', 'ATG', 3)
        b = ReadElement('seq2', 'CCTA', 4)
        print(a + b)
        # seq1seq2 (7): ATGCCTA

        Args:
            other (TypeAlias): ReadElement on the right side of + operator.

        Returns:
            TypeAlias: New ReadElement with concatenated sequences.
        """
        return type(self)(
            self.label + other.label,
            self.sequence + other.sequence,
        )

    def __str__(self) -> str:
        """Print this string when we try to print(object).

        Returns:
            str: Formatted string representation.
        """
        return f'{self.label} ({self.length}): {self.sequence}'

    def __repr__(self) -> str:
        """Print this string when we try to print a list of these objects.

        Returns:
            str: Formatted string representation.
        """
        return self.__str__()

    def __eq__(self, other: TypeAlias = Self) -> bool:
        """Check if two objects are equal.

        Args:
            other (TypeAlias): Object on RHS of equality.

        Returns:
            bool: True if label sequence and length are the same.
        """
        return all((
            self.label == other.label,
            self.sequence == other.sequence,
            self.length == other.length,
        ))


def concat_read_elements(read_elements: list[ReadElement]) -> ReadElement:
    return sum(read_elements, start=ReadElement())
