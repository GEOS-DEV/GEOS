import numpy as np
import hdf5_wrapper


def print_database_iterative(database, level=0):
    """
  @brief print the database targets iteratively by level
  @param database the wrapper for the current database
  @param level the depth within the database
  """
    # Note: you can also iterate over the hdf5_wrapper object directly
    for k in database.keys():
        print("%s%s" % ("  " * level, k))

        if isinstance(database[k], hdf5_wrapper.hdf5_wrapper):
            # This is a group, so continue iterating downward
            print_database_iterative(database[k], level + 1)
        else:
            # This is likely to be an array
            print(database[k])
            print()


def read_write_hdf5_database_example():
    """
  @brief simple demonstration of hdf5_wrapper
  """

    # ------------------------
    # Generate test data
    # ------------------------
    source_a = {
        "1D_double_array": np.random.randn(10),
        "string_array": np.array(["a", "list", "of", "strings"]),
        "child_a": {"2D_double_array": np.random.randn(2, 3)},
    }

    source_b = {
        "1D_integer_array": np.random.randint(0, 100, 5),
        "child_b": {"3D_double_array": np.random.randn(4, 5, 2)},
    }

    # ------------------------
    # Write databases to file
    # ------------------------
    # Write the first piece-by-piece to an hdf5_file
    # Note: when you exit the following scope, the database is automatically closed
    with hdf5_wrapper.hdf5_wrapper("database_a.hdf5", mode="a") as database_a:
        # Assign the two array objects to this level
        database_a["1D_double_array"] = source_a["1D_double_array"]
        database_a["string_array"] = source_a["string_array"]

        # Create a child group and assign the final array
        child_a = database_a["child_a"]
        child_a["2D_double_array"] = source_a["child_a"]["2D_double_array"]

    # Automatically write the second source to a second database
    with hdf5_wrapper.hdf5_wrapper("database_b.hdf5", mode="a") as database_b:
        database_b["/"] = source_b

    # Create a third database that links the either two
    with hdf5_wrapper.hdf5_wrapper("database_c.hdf5", mode="a") as database_c:
        database_c.link("database_a", "database_a.hdf5")
        database_c.link("database_b", "database_b.hdf5")

    # ---------------------------------------
    # Read the databases from the filesystem
    # ---------------------------------------
    print("Database contents:")
    with hdf5_wrapper.hdf5_wrapper("database_c.hdf5") as database_c:
        # Iteratively print the database contents
        print_database_iterative(database_c, 1)

        # As a final note, you can also access low-level h5py functionality
        # by interacting directly with the database target, e.g.:
        print("Database attributes:")
        print("  ", database_c.target.attrs)


if __name__ == "__main__":
    read_write_hdf5_database_example()
