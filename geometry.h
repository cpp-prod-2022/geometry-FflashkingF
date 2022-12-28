#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>
#include <utility>

const double EPS1 = 1e-8;
const double PI = acos(-1);

bool equal(double a, double b) { return std::fabs(a - b) < EPS1; }

bool smaller_or_equal(double a, double b) { return a <= b || equal(a, b); }

bool sign(double a) {
  if (equal(a, 0) || a >= 0) {
    return 1;
  } else {
    return 0;
  }
}

double to_radian(double angle) { return angle * PI / 180; }

namespace {
struct Vector;
}

struct Point {
  double x, y;
  Point() {}
  Point(double x, double y) : x(x), y(y) {}
  explicit Point(const Vector&);
};

std::ostream& operator<<(std::ostream& out, const Point& pt) {
  out << pt.x << ' ' << pt.y << std::endl;
  return out;
}

double dist(const Point& pt1, const Point& pt2) {
  return hypot(pt1.x - pt2.x, pt1.y - pt2.y);
}

bool point_between(const Point& left, const Point& pt, const Point& right) {
  return equal(dist(left, right), dist(left, pt) + dist(pt, right));
}

bool operator==(const Point& pt1, const Point& pt2) {
  return equal(pt1.x, pt2.x) && equal(pt1.y, pt2.y);
}

bool operator!=(const Point& a, const Point& b) { return !(a == b); }

namespace {
struct Vector {
  double x, y;

  Vector() {}
  explicit Vector(double x, double y) : x(x), y(y) {}
  explicit Vector(const Point& start, const Point& end)
      : x(end.x - start.x), y(end.y - start.y) {}
  explicit Vector(const Point& pt) : x(pt.x), y(pt.y) {}

  void rotate(double angle) {
    double cosa = cos(angle);
    double sina = sin(angle);
    double newx = x * cosa - y * sina;
    y = x * sina + y * cosa;
    x = newx;
  }

  void rotate90() {
    std::swap(x, y);
    x *= -1;
  }

  Vector operator-() const { return Vector(-x, -y); }

  Vector& operator*=(double h) {
    x *= h;
    y *= h;
    return *this;
  }

  double len() const { return hypot(x, y); }

  void normalize() {
    double len = this->len();
    x /= len;
    y /= len;
  }
};
}  // namespace

/*std::ostream& operator<<(std::ostream& out, const Vector& v) {
  out << v.x << ' ' << v.y << std::endl;
  return out;
}*/

Point::Point(const Vector& v) : x(v.x), y(v.y) {}

/*Vector operator+(const Vector& v1, const Vector& v2) {
  return Vector(v1.x + v2.x, v1.y + v2.y);
}*/

Vector operator-(const Vector& v1, const Vector& v2) {
  return Vector(v1.x - v2.x, v1.y - v2.y);
}

Point operator+(const Point& pt, const Vector& v) {
  return Point(pt.x + v.x, pt.y + v.y);
}

Point operator-(const Point& pt, const Vector& v) {
  return Point(pt.x - v.x, pt.y - v.y);
}

double operator*(const Vector& v1, const Vector& v2) {
  return v1.x * v2.x + v1.y * v2.y;
}

double operator%(const Vector& v1, const Vector& v2) {
  return v1.x * v2.y - v2.x * v1.y;
}

double angle_between(const Vector& v1, const Vector& v2) {
  return atan2(v1 % v2, v1 * v2);
}

class Line {
 private:
  double a, b, c;

 public:
  Line(const Point& pt0, const Point& pt1)
      : a(pt1.y - pt0.y),
        b(pt0.x - pt1.x),
        c(pt0.y * (pt1.x - pt0.x) + pt0.x * (pt0.y - pt1.y)) {
    // std::cerr << "HOW?" << std::endl;
    // std::cerr << pt0;
    // std::cerr << pt1;
  }
  Line(double k, double b) : a(k), b(-1), c(b) {
    // std::cerr << "NOWAY1" << std::endl;
  }
  Line(const Point& pt, double k) : a(k), b(-1), c(pt.y - k * pt.x) {
    // std::cerr << "NOWAY2" << std::endl;
  }
  //explicit Line(const Point& pt, const Vector& v) : Line(pt, pt + v) {
    // std::cerr << "LINE CONSTR pt v" << std::endl;
    // std::cerr << pt;
    // std::cerr << v;
  //}

  bool operator==(const Line& l) const {
    // if(!equal(a * l.b, b * l.a)) return 0;
    if (!equal(a, 0) || fabs(a) > fabs(b)) {
      double k = l.a / a;
      if (equal(l.b, b * k) && equal(l.c, c * k)) {
        return 1;
      } else {
        return 0;
      }
    } else {
      double k = l.b / b;
      if (equal(l.a, a * k) && equal(l.c, c * k)) {
        return 1;
      } else {
        return 0;
      }
    }
  }

  bool operator!=(const Line& l) const { return !(*this == l); }

  Vector norm_vect() const {
    double len = hypot(a, b);
    return Vector(a / len, b / len);
  }

  double dist(const Point& pt) const {
    return fabs(pt.x * a + pt.y * b + c) / hypot(a, b);
  }

  Point projection(const Point& pt) const {
    Vector v = norm_vect();
    v *= dist(pt);
    Point pt1 = pt + v;
    Point pt2 = pt - v;
    if (dist(pt1) < dist(pt2)) {
      return pt1;
    } else {
      return pt2;
    }
  }

  Point reflect(const Point& pt) const {
    Point projection = this->projection(pt);
    return projection + Vector(pt, projection);
  }

  Point intersection(const Line& l2) const {
    // std::cerr << "intersection" << std::endl;
    // std::cerr << *this;
    // std::cerr << l2;
    double d = a * l2.b - l2.a * b;
    return Point((b * l2.c - l2.b * c) / d, (c * l2.a - a * l2.c) / d);
  }

  friend std::ostream& operator<<(std::ostream&, const Line&);
};

std::ostream& operator<<(std::ostream& out, const Line& l) {
  out << l.a << ' ' << l.b << ' ' << l.c << std::endl;
  return out;
}

class Shape {
 public:
  virtual double perimeter() const = 0;
  virtual double area() const = 0;
  virtual bool containsPoint(const Point&) const = 0;
  virtual bool is_equal(const Shape&) const = 0;
  virtual bool isCongruentTo(const Shape&) const = 0;
  virtual bool isSimilarTo(const Shape&) const = 0;
  virtual void rotate(const Point&, double) = 0;
  virtual void reflect(const Point&) = 0;
  virtual void reflect(const Line&) = 0;
  virtual void scale(const Point&, double) = 0;

  virtual ~Shape() = default;
};

bool operator==(const Shape& sh1, const Shape& sh2) {
  return sh1.is_equal(sh2);
}

bool operator!=(const Shape& sh1, const Shape& sh2) {
  return !sh1.is_equal(sh2);
}

class Ellipse : public Shape {
 protected:
  Point f1, f2;
  double a;

 public:
  Ellipse(const Point& f1, const Point& f2, double _2a)
      : f1(f1), f2(f2), a(_2a / 2) {}

  std::pair<Point, Point> focuses() const { return std::make_pair(f1, f2); }

  double get_c() const { return hypot(f2.x - f1.x, f2.y - f1.y) / 2; }

  double get_b() const {
    double c = get_c();
    return sqrt(a * a - c * c);
  }

  double eccentricity() const { return get_c() / a; }

  Point center() const { return Point((f1.x + f2.x) / 2, (f1.y + f2.y) / 2); }

  std::pair<Line, Line> directrices() const {
    Point center(this->center());
    Vector v(center, f1);

    Vector norm = v;
    std::swap(norm.x, norm.y);
    norm.x *= -1;

    double k = a / get_c();
    v *= k * k;
    return {Line(center + v, center + v + norm), Line(center - v, center - v + norm)};
  }

  double perimeter() const override {  //////////
   // return 1;       
    double b = get_b();        
    //return PI * (3 * (a + b) - sqrt((3 * a + b) * (a + 3 * b)));           ////////////////////////////
    return 4 * a * std::comp_ellint_2(eccentricity());
  }  ///////////////

  double area() const override { return PI * a * get_b(); }

  bool operator==(const Ellipse& e) const {
    return equal(a, e.a) &&
           ((f1 == e.f1 && f2 == e.f2) || (f1 == e.f2 && f2 == e.f1));
  }

  bool operator!=(const Ellipse& e) { return !(*this == e); }

  bool is_equal(const Shape& sh) const final {
    const Ellipse* pointer = dynamic_cast<const Ellipse*>(&sh);
    if (pointer) {
      return *this == *pointer;
    } else {
      return 0;
    }
  }

  bool isCongruentTo(const Shape& sh) const final {
    const Ellipse* e = dynamic_cast<const Ellipse*>(&sh);
    if (e) {
      return a == e->a && equal(get_c(), e->get_c());
    } else {
      return 0;
    }
  }

  bool isSimilarTo(const Shape& sh) const final {
    const Ellipse* e = dynamic_cast<const Ellipse*>(&sh);
    if (e) {
      return equal(a * e->get_c(), get_c() * e->a);
    } else {
      return 0;
    }
  }

  bool containsPoint(const Point& pt) const final {
    if (smaller_or_equal(dist(f1, pt) + dist(f2, pt), 2 * a)) {
      return 1;
    } else {
      return 0;
    }
  }

  void rotate(const Point& center, double angle) final {
    angle = to_radian(angle);
    Vector v1(center, f1), v2(center, f2);
    v1.rotate(angle);
    v2.rotate(angle);
    f1 = center + v1;
    f2 = center + v2;
  }

  void reflect(const Point& center) final {
    f1 = center + Vector(f1, center);
    f2 = center + Vector(f2, center);
  }

  void reflect(const Line& axis) final {
    f1 = axis.reflect(f1);
    f2 = axis.reflect(f2);
  }

  void scale(const Point& pt, double coefficient) final {
    if (f1 == f2) {
      Vector v(pt, f1);
      v *= coefficient;
      f1 = f2 = pt + v;
      a *= coefficient;
    } else {
      Vector v1(pt, f1);
      Vector v2(pt, f2);
      Point center_of_ellipse = this->center();
      Vector v3(center_of_ellipse, f1);
      v3 *= a / v3.len();
      Point pt_on_ellipse = center_of_ellipse + v3;
      Vector v4(pt, pt_on_ellipse);
      v1 *= coefficient;
      v2 *= coefficient;
      v4 *= coefficient;
      f1 = Point(v1);
      f2 = Point(v2);
      a = ((v1 - v4).len() + (v2 - v4).len()) / 2;
    }
  }
};

class Circle : public Ellipse {
 public:
  Circle(const Point& center, double r) : Ellipse(center, center, r * 2) {}
  double radius() const { return a; }

  double area() const final { return PI * a * a; }

  double perimeter() const final { return 2 * PI * a; }
};

class Polygon : public Shape {
 protected:
  mutable std::vector<Point> vert;

 public:
  Polygon(const std::vector<Point>& a) : vert(a) {}
  Polygon() {}
  template <typename... Why>  //??????
  Polygon(Why... args) {      // Point... why not this?
    (vert.emplace_back(args), ...);
  }
  Polygon(const std::initializer_list<Point>& list) {
    for(const Point& i : list) {
      vert.push_back(i);
    }
  }

  size_t verticesCount() const { return vert.size(); }

  const std::vector<Point>& getVertices() const { return vert; }

  bool isConvex() const {
    bool start_sign = sign(Vector(vert.back(), vert[0]) % Vector(vert[0], vert[1]));
    for (size_t i = 1; i < vert.size(); ++i) {
      if (sign(Vector(vert[i - 1], vert[i]) % Vector(vert[i], vert[next(i)])) !=
          start_sign) {
        return 0;
      }
    }
    return 1;
  }

  double perimeter() const final {
    double ans = dist(vert.back(), vert[0]);
    for (size_t i = 1; i < vert.size(); ++i) {
      ans += dist(vert[i], vert[i - 1]);
    }
    return ans;
  }

  double area() const final {
    double ans = 0;
    Vector v1, v2(vert[0], vert[1]);
    for (size_t i = 2; i < vert.size(); ++i) {
      std::swap(v1, v2);
      v2 = Vector(vert[0], vert[i]);
      ans += v1 % v2;
    }
    return fabs(ans) / 2;
  }

  bool help_to_equal(const Polygon& p) const {
    size_t start = static_cast<size_t>(-1);
    for (size_t i = 0; i < vert.size(); ++i) {
      if (vert[i] == p.vert[0]) {
        start = i;
        break;
      }
    }
    if (start != static_cast<size_t>(-1)) {
      std::vector<Point> temp(vert.size());
      copy(vert.begin() + static_cast<int>(start), vert.end(), temp.begin());
      copy(vert.begin(), vert.begin() + static_cast<int>(start),
           temp.begin() + static_cast<int>(vert.size() - start));
      vert = temp;
      for (size_t i = 1; i < vert.size(); ++i) {
        if (vert[i] != p.vert[i]) {
          return 0;
        }
      }
      return 1;
    }
    return 0;
  }

  bool operator==(const Polygon& p) const {
    if (vert.size() != p.vert.size()) return 0;
    bool eq = help_to_equal(p);
    if (eq) return 1;
    std::reverse(vert.begin(), vert.end());
    eq = help_to_equal(p);
    return eq;
  }

  bool operator!=(const Polygon& p) const { return !(*this == p); }

  bool is_equal(const Shape& sh) const final {
    const Polygon* pointer = dynamic_cast<const Polygon*>(&sh);
    if (pointer) {
      return *this == *pointer;
    } else {
      return 0;
    }
  }

  size_t before(size_t i) const { return i == 0 ? vert.size() - 1 : i - 1; }

  size_t next(size_t i) const { return i == vert.size() - 1 ? 0 : i + 1; }

  bool help_to_congruent(const Polygon& p) const {
    for (size_t start = 0; start < vert.size(); ++start) {
      bool ok = 1;
      for (size_t i = start, j = 0; j < vert.size(); ++j, i = next(i)) {
        Vector v1(vert[before(i)], vert[i]), v2(vert[i], vert[next(i)]);
        Vector v3(p.vert[before(j)], p.vert[j]), v4(p.vert[j], p.vert[next(j)]);
        // double a1 = std::fabs(angle_between(v1, v2));
        // double a2 = std::fabs(angle_between(v3, v4));
        if (!equal(v1.len(), v3.len()) || !equal(v2.len(), v4.len()) ||
            !equal(std::fabs(angle_between(v1, v2)),
                   std::fabs(angle_between(v3, v4)))) {
          ok = 0;
          break;
        }
      }
      if (ok) return 1;
    }
    return 0;
  }

  bool isCongruentTo(const Polygon& p) const {
    if (vert.size() != p.vert.size()) return 0;
    bool cong = help_to_congruent(p);
    if (cong) return 1;
    std::reverse(vert.begin(), vert.end());
    cong = help_to_congruent(p);
    return cong;
  }

  bool isCongruentTo(const Shape& sh) const final {
    const Polygon* pointer = dynamic_cast<const Polygon*>(&sh);
    if (pointer) {
      return isCongruentTo(*pointer);
    } else {
      return 0;
    }
  }

  bool help_to_similar(const Polygon& p) const {
    for (size_t start = 0; start < vert.size(); ++start) {
      bool ok = 1;
      double k = -1;
      for (size_t i = start, j = 0; j < vert.size(); ++j, i = next(i)) {
        // std::cerr << "i: " << i << std::endl;
        // std::cerr << "j: " << j << std::endl;
        Vector v1(vert[before(i)], vert[i]), v2(vert[i], vert[next(i)]);
        Vector v3(p.vert[before(j)], p.vert[j]), v4(p.vert[j], p.vert[next(j)]);
        if (j == 0) {
          k = v1.len() / v3.len();
          // std::cerr << "k: " << k << std::endl;
        }
        v3 *= k;
        v4 *= k;
        double a1 = std::fabs(angle_between(v1, v2));
        double a2 = std::fabs(angle_between(v3, v4));

        // std::cerr << "vec1: " << v1.len() << std::endl;
        // std::cerr << "vec2: " << v2.len() << std::endl;
        // std::cerr << "vec3: " << v3.len() << std::endl;
        // std::cerr << "vec4: " << v4.len() << std::endl;
        // std::cerr << "a1: " << a1 << std::endl;
        // std::cerr << "a2: " << a2 << std::endl;
        if (!equal(v1.len(), v3.len()) || !equal(v2.len(), v4.len()) ||
            !equal(a1, a2)) {
          ok = 0;
          // std::cerr << "shish" << std::endl;
          break;
        }
      }
      if (ok) return 1;
    }
    return 0;
  }

  bool isSimilarTo(const Polygon& p) const {
    // std::cerr << "isSim" << std::endl;
    // std::cerr << "OFSIJFOSDJIF" << std::endl;
    if (vert.size() != p.vert.size()) return 0;
    // std::cerr << "sz == p.sz" << std::endl;
    bool sim = help_to_similar(p);
    if (sim) return 1;
    std::reverse(vert.begin(), vert.end());
    sim = help_to_similar(p);
    // if(sim)//std::cerr << "YES2" << std::endl;
    // else//std::cerr << "NO2" << std::endl;
    return sim;
  }

  bool isSimilarTo(const Shape& sh) const final {
    const Polygon* pointer = dynamic_cast<const Polygon*>(&sh);
    if (pointer) {
      // std::cerr << "BEN" << std::endl;
      return isSimilarTo(*pointer);
    } else {
      return 0;
    }
  }

  bool containsPoint(const Point& pt) const final {
    if (point_between(vert[0], pt, vert.back())) {
      return 1;
    }
    for (size_t i = 1; i < vert.size(); ++i) {
      if (point_between(vert[i - 1], pt, vert[i])) {
        return 1;
      }
    }
    double param = angle_between(Vector(pt, vert.back()), Vector(pt, vert[0]));
    Vector v1, v2(pt, vert[0]);
    for (size_t i = 1; i < vert.size(); ++i) {
      std::swap(v1, v2);
      v2 = Vector(pt, vert[i]);
      param += angle_between(v1, v2);
    }
    return !equal(param, 0);
  }

  void rotate(const Point& center, double angle) final {
    angle = to_radian(angle);
    for (size_t i = 0; i < vert.size(); ++i) {
      Vector v(center, vert[i]);
      v.rotate(angle);
      vert[i] = center + v;
    }
  }

  void reflect(const Point& center) final {
    for (size_t i = 0; i < vert.size(); ++i) {
      Vector v(vert[i], center);
      vert[i] = center + Vector(vert[i], center);
      ;
    }
  }

  void reflect(const Line& axis) final {
    for (size_t i = 0; i < vert.size(); ++i) {
      vert[i] = axis.reflect(vert[i]);
    }
  }

  void scale(const Point& center, double coefficeint) final {
    for (size_t i = 0; i < vert.size(); ++i) {
      Vector v(center, vert[i]);
      v *= coefficeint;
      vert[i] = center + v;
    }
  }

  friend std::ostream& operator<<(std::ostream&, const Polygon&);
};

std::ostream& operator<<(std::ostream& out, const Polygon& p) {
  for (size_t i = 0; i < p.vert.size(); ++i) {
    out << p.vert[i].x << ' ' << p.vert[i].y << std::endl;
  }
  return out;
}

class Rectangle : public Polygon {
 public:
  Rectangle(const Point& pt1, const Point& pt2, double k) {
    if (k < 1) {
      k = 1 / k;
    }
    Vector v(pt1, pt2);
    double coeff = 1 / hypot(1, k);
    double angle = acos(coeff);
    Vector v2 = v;
    v2.rotate(angle);
    v2 *= coeff;
    vert = {pt1, pt1 + v2, pt2, pt2 - v2};
  }

  Point center() const {
    return Point(vert[0].x + vert[2].x, vert[0].y + vert[2].y);
  }

  std::pair<Line, Line> diagonals() const {
    return {Line(vert[0], vert[2]), Line(vert[1], vert[3])};
  }
};

class Square : public Rectangle {
 public:
  Square(const Point& pt1, const Point& pt2) : Rectangle(pt1, pt2, 1) {}

  Circle inscribedCircle() const {
    return Circle(center(), dist(vert[0], vert[1]) / 2);
  }

  Circle circumscribedCircle() const {
    return Circle(center(), dist(vert[0], vert[2]) / 2);
  }
};

class Triangle : public Polygon {
 public:
  Triangle(const Point& a, const Point& b, const Point& c) : Polygon(a, b, c) {}

  Point centroid() const {
    return Point((vert[0].x + vert[1].x + vert[2].x) / 3,
                 (vert[0].y + vert[1].y + vert[2].y) / 3);
  }

  Point orthocenter() const {
    Line l1(vert[1], vert[2]);
    Line h1(vert[0], l1.projection(vert[0]));
    Line l2(vert[2], vert[0]);
    Line h2(vert[1], l2.projection(vert[1]));
    return h1.intersection(h2);
  }

  Circle circumscribedCircle() const {
    // std::cerr << "cricum" << std::endl;
    Point c1((vert[0].x + vert[1].x) / 2, (vert[0].y + vert[1].y) / 2);
    Vector v1(vert[0], vert[1]);
    v1.rotate90();
    Line sp1(c1, c1 + v1);
    // std::cerr << "c1: " << c1;
    // std::cerr << "v1: " << v1;
    // std::cerr << "Line sp1: " << sp1;

    Point c2((vert[1].x + vert[2].x) / 2, (vert[1].y + vert[2].y) / 2);
    Vector v2(vert[1], vert[2]);
    v2.rotate90();
    Line sp2(c2, c2 + v2);
    // std::cerr << "Line sp2: " << sp2;

    Point center = sp1.intersection(sp2);
    // std::cerr << "ceter: " << center;
    // std::cerr << "vert[0]: " << vert[0];
    // std::cerr << "radius: " << dist(center, vert[0]) << std::endl;
    return Circle(center, dist(center, vert[0]));
  }

  Circle inscribedCircle() const {
    // std::cerr << "inscribedCircle" << std::endl;
    Vector v1(vert[0], vert[1]);
    Vector v2(vert[0], vert[2]);
    v1.normalize();
    v2.normalize();
    // std::cerr << "first" << std::endl;
    // std::cerr << vert[0];
    // std::cerr << v1 + v2;
    Line l1(vert[0], vert[0] + v1 + v2);
    // std::cerr << "line1: " << l1;

    v1 = Vector(vert[1], vert[2]);
    v2 = Vector(vert[1], vert[0]);
    v1.normalize();
    v2.normalize();
    // std::cerr << "second" << std::endl;
    // std::cerr << vert[1];
    // std::cerr << v1 + v2;
    Line l2(vert[1], vert[1] + v1 + v2);

    Point center(l1.intersection(l2));
    // std::cerr << "area: " << area() << std::endl;
    // std::cerr << "perimeter: " << perimeter() << std::endl;
    // std::cerr << std::endl;
    return Circle(center, area() * 2 / perimeter());
  }

  Line EulerLine() const { return Line(centroid(), orthocenter()); }

  Circle ninePointsCircle() const {
    // std::cerr << std::endl << std::endl <<  "ninePointsCirce" << std::endl;
    Point pt1((vert[0].x + vert[1].x) / 2, (vert[0].y + vert[1].y) / 2);
    Point pt2((vert[0].x + vert[2].x) / 2, (vert[0].y + vert[2].y) / 2);
    Point pt3((vert[2].x + vert[1].x) / 2, (vert[2].y + vert[1].y) / 2);

    Vector v1(pt1, pt2);
    Point c1((pt1.x + pt2.x) / 2, (pt1.y + pt2.y) / 2);
    v1.rotate90();
    Line l1(c1, c1 + v1);

    Vector v2(pt1, pt3);
    Point c2((pt1.x + pt3.x) / 2, (pt1.y + pt3.y) / 2);
    v2.rotate90();
    Line l2(c2, c2 + v2);

    Point center(l1.intersection(l2));
    return Circle(center, dist(center, pt1));
  }

  // friend std::ostream& operator<<(std::ostream&, const Circle&);
};

std::ostream& operator<<(std::ostream& out, const Circle& c) {
  out << "start" << std::endl;
  out << c.focuses().first << c.radius() << std::endl;
  out << "endl" << std::endl;
  return out;
}
