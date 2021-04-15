class Trajectories:
    ONE_MAX = 1
    FIRST_HIGHER = 2
    SECOND_HIGHER = 3

    @staticmethod
    def classify(maxima):
        if maxima.shape[0] < 2:
            return Trajectories.ONE_MAX
        if maxima.iloc[0]["dim"] > maxima.iloc[1]["dim"]:
            return Trajectories.FIRST_HIGHER
        return Trajectories.SECOND_HIGHER

    @staticmethod
    def iter():
        return (
            (Trajectories.ONE_MAX, {"color": "tab:green", "mean_color": "tab:brown"}),
            (
                Trajectories.FIRST_HIGHER,
                {"color": "tab:blue", "mean_color": "tab:purple"},
            ),
            (
                Trajectories.SECOND_HIGHER,
                {"color": "tab:orange", "mean_color": "tab:red"},
            ),
        )
