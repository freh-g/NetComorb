// This file is a slight modification of page_rank.hpp in Boost Graph Library
// to calculate Page Rank with weighted priors. Thus original license of the 
// code is inhereted for the code (see below for original license).
// 2011 - Emre Guney (Universitat Pompeu Fabra)
//
// Copyright 2004-5 The Trustees of Indiana University.
// Copyright 2002 Brad King and Douglas Gregor

// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

//  Authors: Douglas Gregor
//           Andrew Lumsdaine

#ifndef BOOST_GRAPH_PAGE_RANK_WITH_WEIGHTED_PRIORS_HPP
#define BOOST_GRAPH_PAGE_RANK_WITH_WEIGHTED_PRIORS_HPP

#include <boost/property_map/property_map.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/properties.hpp>
#include <boost/graph/iteration_macros.hpp>
#include <vector>

//#include <iostream>

namespace boost { namespace graph {

struct n_iterations
{
  explicit n_iterations(std::size_t n) : n(n) { }

  template<typename RankMap, typename Graph>
  bool 
  operator()(const RankMap&, const Graph&)
  {
    return n-- == 0;
  }

 private:
  std::size_t n;
};

namespace detail {
  template<typename Graph, typename RankMap, typename RankMap2, typename RankMap3, typename WeightMap>
  void page_rank_with_weighted_priors_step(const Graph& g, RankMap from_rank, RankMap2 to_rank, RankMap3 initial_rank_map, WeightMap weight,
                      typename property_traits<RankMap>::value_type damping,
		      typename graph_traits<Graph>::vertices_size_type n,
                      incidence_graph_tag)
  {
    typedef typename property_traits<RankMap>::value_type rank_type;

    // Set new rank maps 
    //BGL_FORALL_VERTICES_T(v, g, Graph) put(to_rank, v, rank_type((1 - damping)/n));
    BGL_FORALL_VERTICES_T(v, g, Graph) put(to_rank, v, rank_type((1 - damping) * get(initial_rank_map, v)));
 
    BGL_FORALL_VERTICES_T(u, g, Graph) {
      //rank_type u_rank_out = damping * get(from_rank, u) / out_degree(u, g);
      rank_type u_rank_out = damping * get(from_rank, u) / out_degree(u, g);
      BGL_FORALL_ADJ_T(u, v, g, Graph)
        put(to_rank, v, get(to_rank, v) + get(weight, edge(u, v, g)) * u_rank_out);
    }
  }

  template<typename Graph, typename RankMap, typename RankMap2, typename RankMap3, typename WeightMap>
  void page_rank_with_weighted_priors_step(const Graph& g, RankMap from_rank, RankMap2 to_rank, RankMap3 initial_rank_map, WeightMap weight,
                      typename property_traits<RankMap>::value_type damping, 
		      typename graph_traits<Graph>::vertices_size_type n,
                      bidirectional_graph_tag)
  {
    typedef typename property_traits<RankMap>::value_type damping_type;

    BGL_FORALL_VERTICES_T(v, g, Graph) {
      typename property_traits<RankMap>::value_type rank(0);
      //BGL_FORALL_INEDGES_T(v, e, g, Graph)
      //  rank += get(from_rank, source(e, g)) / out_degree(source(e, g), g);
      BGL_FORALL_INEDGES_T(v, e, g, Graph)
        rank += get(weight, e) * get(from_rank, source(e, g)) / out_degree(source(e, g), g);
	//std::cout << "accumulating rank " << v << ": " << rank << std::endl;
	//std::cout << "final value: " << (damping_type(1) - damping) / n + damping * rank << std::endl;
      //put(to_rank, v, (damping_type(1) - damping) / n + damping * rank);
      put(to_rank, v, (damping_type(1) - damping) * get(initial_rank_map, v) + damping * rank);
    }
  }
} // end namespace detail

template<typename Graph, typename RankMap, typename WeightMap, typename Done, typename RankMap2>
void
page_rank_with_weighted_priors(const Graph& g, RankMap initial_rank_map, WeightMap weight, RankMap rank_map, 
          Done done, 
          typename property_traits<RankMap>::value_type damping,
          typename graph_traits<Graph>::vertices_size_type n,
          RankMap2 rank_map2)
{
  typedef typename property_traits<RankMap>::value_type rank_type;

  //rank_type initial_rank = rank_type(rank_type(1) / n);
  //BGL_FORALL_VERTICES_T(v, g, Graph) put(rank_map, v, initial_rank);
  // Assigning given initial ranks at the begining of the algorithm as pageranks at step 0 - has insignificant effect
  BGL_FORALL_VERTICES_T(v, g, Graph) put(rank_map, v, rank_type(get(initial_rank_map, v)));


  bool to_map_2 = true;
  while ((to_map_2 && !done(rank_map, g)) ||
         (!to_map_2 && !done(rank_map2, g))) {
    typedef typename graph_traits<Graph>::traversal_category category;

    if (to_map_2) {
      detail::page_rank_with_weighted_priors_step(g, rank_map, rank_map2, initial_rank_map, weight, damping, n, category());
    } else {
      detail::page_rank_with_weighted_priors_step(g, rank_map2, rank_map, initial_rank_map, weight, damping, n, category());
    }
    to_map_2 = !to_map_2;
  }

  if (!to_map_2) {
    BGL_FORALL_VERTICES_T(v, g, Graph) put(rank_map, v, get(rank_map2, v));
  }
}

template<typename Graph, typename RankMap, typename WeightMap, typename Done>
void
page_rank_with_weighted_priors(const Graph& g, RankMap initial_rank_map, WeightMap weight, RankMap rank_map, 
          Done done, 
          typename property_traits<RankMap>::value_type damping,
          typename graph_traits<Graph>::vertices_size_type n)
{
  typedef typename property_traits<RankMap>::value_type rank_type;

  std::vector<rank_type> ranks2(num_vertices(g));
  page_rank_with_weighted_priors(g, initial_rank_map, weight, rank_map, done, damping, n,
            make_iterator_property_map(ranks2.begin(), get(vertex_index, g)));
}

template<typename Graph, typename RankMap, typename WeightMap, typename Done>
inline void
page_rank_with_weighted_priors(const Graph& g, RankMap initial_rank_map, WeightMap weight, RankMap rank_map, Done done, 
          typename property_traits<RankMap>::value_type damping = 0.85)
{
  page_rank_with_weighted_priors(g, initial_rank_map, weight, rank_map, done, damping, num_vertices(g));
}

template<typename Graph, typename RankMap, typename WeightMap>
inline void
page_rank_with_weighted_priors(const Graph& g, RankMap initial_rank_map, WeightMap weight, RankMap rank_map)
{
  page_rank_with_weighted_priors(g, initial_rank_map, weight, rank_map, n_iterations(20));
}

// TBD: this could be _much_ more efficient, using a queue to store
// the vertices that should be reprocessed and keeping track of which
// vertices are in the queue with a property map. Baah, this only
// applies when we have a bidirectional graph.
template<typename MutableGraph>
void
remove_dangling_links(MutableGraph& g)
{
  typename graph_traits<MutableGraph>::vertices_size_type old_n;
  do {
    old_n = num_vertices(g);

    typename graph_traits<MutableGraph>::vertex_iterator vi, vi_end;
    for (tie(vi, vi_end) = vertices(g); vi != vi_end; /* in loop */) {
      typename graph_traits<MutableGraph>::vertex_descriptor v = *vi++;
      if (out_degree(v, g) == 0) {
        clear_vertex(v, g);
        remove_vertex(v, g);
      }
    }
  } while (num_vertices(g) < old_n);
}

} } // end namespace boost::graph

#endif // BOOST_GRAPH_PAGE_RANK_WITH_WEIGHTED_PRIORS_HPP
