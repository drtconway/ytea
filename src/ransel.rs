
pub trait Ransel
where
{
    type Rank;
    type Position;

    fn size(&self) -> Self::Position;

    fn count(&self) -> Self::Rank;

    fn rank(&self, x: Self::Position) -> Self::Rank;

    fn rank2(&self, x1: Self::Position, x2: Self::Position) -> (Self::Rank, Self::Rank) {
        (self.rank(x1), self.rank(x2))
    }

    fn select(&self, i: Self::Rank) -> Self::Position;
}
