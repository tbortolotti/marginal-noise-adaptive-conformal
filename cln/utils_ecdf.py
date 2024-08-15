import numpy as np

class EmpiricalCDF2DGrid:
    def __init__(self, X, Y):
        """
        Initialize the 2D empirical CDF by sorting the (X, Y) pairs together.

        Parameters:
        X (np.array): 1D array of observed values for the first random variable
        Y (np.array): 1D array of observed values for the second random variable
        """
        self.n = len(X)
        # Combine X and Y into pairs and sort by X first, then Y (lexicographically)
        sorted_indices = np.lexsort((Y, X))
        self.X_sorted = X[sorted_indices]
        self.Y_sorted = Y[sorted_indices]
        self.grid_evaluated = False

    def evaluate_grid(self, x_grid, y_grid):
        """
        Evaluate the empirical CDF at each point on the provided grid of (x, y) values.

        Parameters:
        x_grid (np.array): 1D array of x-coordinates where the CDF is evaluated
        y_grid (np.array): 1D array of y-coordinates where the CDF is evaluated

        Returns:
        np.array: A 2D array of ECDF values for each (x, y) combination on the grid
        """
        # Store the grid for later querying
        self.x_grid = x_grid
        self.y_grid = y_grid

        # Precompute positions of x_grid and y_grid in the sorted data
        x_positions = np.searchsorted(self.X_sorted, x_grid, side='right')
        cdf_grid = np.zeros((len(x_grid), len(y_grid)))

        for i, x_pos in enumerate(x_positions):
            if x_pos == 0:
                cdf_grid[i, :] = 0.0  # If no points have X <= x_grid[i], the ECDF is 0
            else:
                # Precompute Y counts up to x_pos for the entire y_grid
                y_counts = np.searchsorted(np.sort(self.Y_sorted[:x_pos]), y_grid, side='right')
                cdf_grid[i, :] = y_counts / self.n

        # Mark that the grid has been evaluated
        self.grid_evaluated = True
        self.cdf_grid = cdf_grid

        return cdf_grid

    def query(self, x, y):
        """
        Query the empirical CDF at point (x, y). If the point is on the precomputed grid,
        return the exact grid result; otherwise, calculate using the search method.

        Parameters:
        x (float): The x-coordinate where the CDF is evaluated
        y (float): The y-coordinate where the CDF is evaluated

        Returns:
        float: The value of the empirical CDF at point (x, y)
        """
        # Check if the grid has been evaluated
        if self.grid_evaluated:
            # Try to find the exact match in the grid
            x_idx = np.searchsorted(self.x_grid, x, side='left')
            y_idx = np.searchsorted(self.y_grid, y, side='left')

            # If the point belongs to the grid, return the grid value
            if (x_idx < len(self.x_grid) and y_idx < len(self.y_grid) and
                np.isclose(self.x_grid[x_idx], x) and np.isclose(self.y_grid[y_idx], y)):
                return self.cdf_grid[x_idx, y_idx]

        # If the point is not in the grid, compute it directly using search
        idx_x = np.searchsorted(self.X_sorted, x, side='right')
        if idx_x == 0:
            return 0.0
        count_y = np.sum(self.Y_sorted[:idx_x] <= y)
        return count_y / self.n


# # Example Usage

# # Initialize the 2D empirical CDF with the observed data
# ecdf_2d_grid = EmpiricalCDF2DGrid(X, Y)

# # Define the grid where you want to evaluate the CDF
# x_grid = np.linspace(0, 1, 100)  # 100 points between -3 and 3 for x
# y_grid = np.linspace(0, 1, 100)  # 100 points between -3 and 3 for y

# # Evaluate the ECDF on the grid
# results_grid = ecdf_2d_grid.evaluate_grid(x_grid, y_grid)

# # Test the query function for a grid point and a non-grid point
# grid_query_result = ecdf_2d_grid.query(0.0, 0.0)  # On grid
# non_grid_query_result = ecdf_2d_grid.query(0.15, 0.15)  # Off grid

# # Print the result
# print(f"Grid query result: {grid_query_result}")
# print(f"Non-grid query result: {non_grid_query_result}")
