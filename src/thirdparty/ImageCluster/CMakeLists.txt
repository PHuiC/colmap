COLMAP_ADD_LIBRARY(
    NAME colmap_image_cluster
    SRCS
        lemon_util.h lemon_util.cpp
        graph_types.h
        geo_covi_cluster.h geo_covi_cluster.cpp
        hierarchical_cluster.h hierarchical_cluster.cpp
    PRIVATE_LINK_LIBS
        lemon
		metis
)
