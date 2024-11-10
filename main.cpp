#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Delaunay_mesher_2.h>
#include <CGAL/Delaunay_mesh_face_base_2.h>
#include <CGAL/boost/graph/graph_traits_Constrained_triangulation_2.h>
#include <CGAL/Surface_mesh.h>

#include <iostream>
#include <unordered_map>

#define cimg_display 0
#include "CImg.h"
using namespace cimg_library;

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Triangulation_vertex_base_2<K> Vb;
typedef CGAL::Delaunay_mesh_face_base_2<K> Fb;
typedef CGAL::Triangulation_data_structure_2<Vb, Fb> Tds;
typedef CGAL::Constrained_Delaunay_triangulation_2<K, Tds> CDT;
typedef CGAL::Surface_mesh<K::Point_3> Mesh;

typedef CDT::Vertex_handle Vertex_handle;
typedef CDT::Point Point;


struct Vpm
{
  using key_type = K::Point_2;
  using value_type = K::Point_3;
  using reference_type = value_type;
  using category = boost::readable_property_map_tag;

  friend value_type get(Vpm, const K::Point_2& p2)
  {
    return K::Point_3(p2.x(), p2.y(), 0.);
  }
};


class Criteria :
    public virtual CGAL::Delaunay_mesh_criteria_2<CDT>
{
protected:
  typedef typename CDT::Geom_traits Geom_traits;
  double sizebound;
  CImg<int> image;

public:
  typedef Delaunay_mesh_criteria_2<CDT> Base;

  Criteria(const CImg<int>& image,
           const double aspect_bound = 0.125,
           const double size_bound = 0,
           const Geom_traits& traits = Geom_traits())
    : Base(aspect_bound, traits), sizebound(size_bound)
    , image(image)
    {}

  inline
  double size_bound() const { return sizebound; }

  inline
  void set_size_bound(const double sb) { sizebound = sb; }

  // first: squared_minimum_sine
  // second: size
  struct Quality : public std::pair<double, double>
  {
    typedef std::pair<double, double> Base;

    Quality() : Base() {};
    Quality(double _sine, double _size) : Base(_sine, _size) {}

    const double& size() const { return second; }
    const double& sine() const { return first; }

    // q1<q2 means q1 is prioritized over q2
    // ( q1 == *this, q2 == q )
    bool operator<(const Quality& q) const
    {
      if( size() > 1 )
        if( q.size() > 1 )
          return ( size() > q.size() );
        else
          return true; // *this is big but not q
      else
        if( q.size() >  1 )
          return false; // q is big but not *this
      return( sine() < q.sine() );
    }

    std::ostream& operator<<(std::ostream& out) const
    {
      return out << "(size=" << size()
                 << ", sine=" << sine() << ")";
    }
  };

  class Is_bad: public Base::Is_bad
  {
  protected:
    const double squared_size_bound; // squared size bound on edge length
    const CImg<int>& image;
  public:
    typedef typename Base::Is_bad::Point_2 Point_2;

    Is_bad(const CImg<int>& image,
           const double aspect_bound,
           const double size_bound,
           const Geom_traits& traits)
      : Base::Is_bad(aspect_bound, traits)
      , squared_size_bound(size_bound * size_bound)
      , image(image)
    {}

    CGAL::Mesh_2::Face_badness operator()(const Quality q) const
    {
      if( q.size() > 1 )
        return CGAL::Mesh_2::IMPERATIVELY_BAD;
      if( q.sine() < this->B )
        return CGAL::Mesh_2::BAD;
      else
        return CGAL::Mesh_2::NOT_BAD;
    }

    CGAL::Mesh_2::Face_badness operator()(const typename CDT::Face_handle& fh,
                                    Quality& q) const
    {
      typedef typename CDT::Geom_traits Geom_traits;
      typedef typename Geom_traits::Compute_area_2 Compute_area_2;
      typedef typename Geom_traits::Compute_squared_distance_2
        Compute_squared_distance_2;

      Compute_squared_distance_2 squared_distance =
        this->traits.compute_squared_distance_2_object();

      const Point_2& pa = fh->vertex(0)->point();
      const Point_2& pb = fh->vertex(1)->point();
      const Point_2& pc = fh->vertex(2)->point();

      double
        a = CGAL::to_double(squared_distance(pb, pc)),
        b = CGAL::to_double(squared_distance(pc, pa)),
        c = CGAL::to_double(squared_distance(pa, pb));

      double max_sq_length; // squared max edge length
      double second_max_sq_length;

      if(a<b)
        {
          if(b<c) {
            max_sq_length = c;
            second_max_sq_length = b;
          }
          else { // c<=b
            max_sq_length = b;
            second_max_sq_length = ( a < c ? c : a );
          }
        }
      else // b<=a
        {
          if(a<c) {
            max_sq_length = c;
            second_max_sq_length = a;
          }
          else { // c<=a
            max_sq_length = a;
            second_max_sq_length = ( b < c ? c : b );
          }
        }

      q.second = 0;
      if( squared_size_bound != 0 )
      {
          //          std::cerr << squared_size_bound << std::endl;

          Point_2 cc = CGAL::centroid(fh->vertex(0)->point(),
                                      fh->vertex(1)->point(),
                                      fh->vertex(2)->point());

          int v = image(floor(cc.x()), floor(cc.y()));

          q.second = max_sq_length / squared_size_bound * (0.01+(255-v)/255.);
            // normalized by size bound to deal
            // with size field
          if( q.size() > 1 )
            {
              q.first = 1; // (do not compute sine)
              return CGAL::Mesh_2::IMPERATIVELY_BAD;
            }

      }

      Compute_area_2 area_2 = this->traits.compute_area_2_object();

      double area = 2*CGAL::to_double(area_2(pa, pb, pc));

      q.first = (area * area) / (max_sq_length * second_max_sq_length); // (sine)

      if( q.sine() < this->B )
        return CGAL::Mesh_2::BAD;
      else
        return CGAL::Mesh_2::NOT_BAD;
    }
  };

  Is_bad is_bad_object() const
  { return Is_bad(image, this->bound(), size_bound(),
                  this->traits /* from the bad class */); }
};



int main(int, char** argv)
{
  CImg<int> image(argv[1]);
  std::cout << image.width() << " " << image.height() << "\n";
  int w = image.width();
  int h = image.height();

  std::cout << image(309,309) << "\n";
  std::cout << image(0,0) << "\n";

  CDT cdt;

  Vertex_handle va = cdt.insert(Point(0,0));
  Vertex_handle vb = cdt.insert(Point(0,h));
  Vertex_handle vc = cdt.insert(Point(w,h));
  Vertex_handle vd = cdt.insert(Point(w,0));

  cdt.insert_constraint(va, vb);
  cdt.insert_constraint(vb, vc);
  cdt.insert_constraint(vc, vd);
  cdt.insert_constraint(vd, va);

  CGAL::refine_Delaunay_mesh_2(cdt, CGAL::parameters::criteria(Criteria(image, 0.125, 5)));


  // Mesh sm;
  // CGAL::copy_face_graph(cdt, sm, CGAL::parameters::vertex_point_map(Vpm())); // does not work. Why???
  // std::ofstream("out.off") << sm;

  unsigned int i=0;
  std::unordered_map<CDT::Vertex_handle, unsigned int> v2i;

  std::ofstream out("out.off");
  out << "OFF\n" << cdt.number_of_vertices() << " " << cdt.number_of_faces() << " 0\n";
  for (auto vh : cdt.finite_vertex_handles())
  {
    out << vh->point() << " 0\n";
    v2i[vh]=i++;
  }
  for (auto fh : cdt.finite_face_handles())
  {
    out << "3 " << v2i[fh->vertex(0)] << " " << v2i[fh->vertex(1)] << " " << v2i[fh->vertex(2)] << "\n";
  }

  return 0;
}
