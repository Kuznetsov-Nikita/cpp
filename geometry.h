#include <cmath>
#include <iostream>
#include <vector>
#include <utility>

const double EPS = 1e-4;

struct Point {
public:
    double x, y;

    Point() = default;

    Point(double x, double y): x(x), y(y) {}

    Point& operator+=(const Point another) {
        x += another.x;
        y += another.y;

        return *this;
    }

    Point& operator-=(const Point another) {
        x -= another.x;
        y -= another.y;

        return *this;
    }

    void reflex(const Point point) {
        x = 2.0 * point.x - x;
        y = 2.0 * point.y - y;
    }
};

double operator*(const Point point_1, const Point point_2) {
    return point_1.x * point_2.x + point_1.y * point_2.y;
}

Point operator+(const Point point_1, const Point point_2) {
    Point copy = point_1;
    copy += point_2;

    return copy;
}

Point operator-(const Point point_1, const Point point_2) {
    Point copy = point_1;
    copy -= point_2;

    return copy;
}

double lenght(const Point point_1, const Point point_2) {
    return (sqrt((point_1.x - point_2.x) * (point_1.x - point_2.x) + (point_1.y - point_2.y) * (point_1.y - point_2.y)));
}

bool operator==(const Point& point_1, const Point& point_2) {
    return (fabs(point_1.x - point_2.x) < EPS) && (fabs(point_1.y - point_2.y) < EPS);
}

bool operator!=(const Point& point_1, const Point& point_2) {
    return !(point_1 == point_2);
}

Point get_vec_from_2_points(const Point& point_1, const Point& point_2) {
    return Point ((point_2.x - point_1.x) / lenght(point_1, point_2), (point_2.y - point_1.y) / lenght(point_1, point_2));
}

Point middle_point(const Point& point_1, const Point& point_2) {
    return Point((point_1.x + point_2.x) / 2.0, (point_1.y + point_2.y) / 2.0);
}

double product(const Point& point_1, const Point& point_2, const Point& point_3) {
    return (point_2.x - point_1.x) * (point_3.y - point_1.y) - (point_2.y - point_1.y) * (point_3.x - point_1.x);
}

class Line {
public:
    double a, b, c;

    Line() = default;

    Line(Point point_1, Point point_2) {
        if (fabs(point_1.x - point_2.x) < EPS) {
            a = 0.0;
            b = 1.0;
            c = -point_1.x;
        } else {
            a = 1.0;
            b = -(point_2.y - point_1.y) / (point_2.x - point_1.x);
            c = -(point_1.y + b * point_1.x);
        }
    }

    Line(double k, double b) {
        a = 1.0;
        this->b = -k;
        c = -b;
    }

    Line(Point point, double k) {
        a = 1.0;
        b = -k;
        c = -(point.y + b * point.x);
    }
};

bool operator==(const Line& line_1, const Line& line_2) {
    if (fabs(line_1.a - line_2.a) < EPS && fabs(line_1.b - line_2.b) < EPS && fabs(line_1.c - line_2.c) < EPS) {
        return true;
    } else {
        return false;
    }
}

bool operator!=(const Line& line_1, const Line& line_2) {
    return !(line_1 == line_2);
}


Point l_reflex(Line line, Point point) {
    Point t_point_1, t_point_2;

    if (fabs(line.a) < EPS) {
        t_point_1 = Point(-(line.c / line.b), 0.0);
        t_point_2 = Point(-(line.c / line.b), 1.0);
    } else {
        t_point_1 = Point(0.0, -(line.c / line.a));
        t_point_2 = Point(1.0, -(line.b + line.c) / line.a);
    }

    Point vec = t_point_2 - t_point_1;
    vec = Point(-vec.y / lenght(vec, Point(0.0, 0.0)), vec.x / lenght(vec, Point(0.0, 0.0)));

    double distance = fabs(line.a * point.y + line.b * point.x + line.c) / sqrt(line.a * line.a + line.b * line.b);

    if (line.a * point.y + line.b * point.x + line.c < - EPS) {
        return Point(point.x + vec.x * 2 * distance, point.y + vec.y * 2 * distance);
    } else {
        return Point(point.x - vec.x * 2 * distance, point.y - vec.y * 2 * distance);
    }
}

class Shape {
public:
    virtual double perimeter() const = 0;
    virtual double area() const = 0;

    virtual bool isCongruentTo(const Shape& another) const = 0;
    virtual bool isSimilarTo(const Shape& another) const = 0;
    virtual bool containsPoint(Point point) const = 0;
    virtual void rotate(Point center, double angle) = 0;
    virtual void reflex(Point center) = 0;
    virtual void reflex(Line axis) = 0;
    virtual void scale(Point center, double coefficient) = 0;
    virtual bool operator==(const Shape& another) const = 0;
    virtual bool operator!=(const Shape& another) const = 0;

    virtual ~Shape() = 0;
};

Shape::~Shape() {}

class Ellipse: public Shape {
protected:
    std::pair<Point, Point> f;
    double a = 0.0;

public:
    Ellipse() = default;

    Ellipse(const Point focus_1, const Point focus_2, double distance_sum) {
        f = std::pair<Point, Point> (focus_1, focus_2);
        a = distance_sum / 2.0;
    }

    std::pair<Point, Point> focuses() const {
        return f;
    }

    double get_a() const {
        return a;
    }

    double eccentricity() const {
        return lenght(f.first, f.second) / (2.0 * a);
    }

    Point center() const {
        return middle_point(f.first, f.second);
    }

    std::pair<Line, Line> directrices() const {
        Point vec1 = get_vec_from_2_points(f.second, f.first);
        Point vec2 = Point(-vec1.y, vec1.x);

        Point d1_beg (this->center().x - a / this->eccentricity() * vec1.x, this->center().y - a / this->eccentricity() * vec1.y);
        Point d2_beg (this->center().x + a / this->eccentricity() * vec1.x, this->center().y + a / this->eccentricity() * vec1.y);

        return std::pair<Line, Line> (Line(d1_beg, d1_beg + vec2), Line(d2_beg, d2_beg + vec2));
    }

    double perimeter() const override {
        double b = a * sqrt(1 - lenght(f.first, f.second) / (2.0 * a) * lenght(f.first, f.second) / (2.0 * a));

        return M_PI * (a + b) * (1 + (3 * (a - b) / (a + b) * (a - b) / (a + b)) / (10 + sqrt(4 - 3 * (a - b) / (a + b) * (a - b) / (a + b))));
    }

    double area() const override {
        if (fabs(a) < EPS) {
            return 0.0;
        }

        double b = a * sqrt(1 - lenght(f.first, f.second) / (2.0 * a) * lenght(f.first, f.second) / (2.0 * a));
        return M_PI * a * b;
    }

    virtual bool isCongruentTo(const Shape& another) const override {
        try {
            Ellipse another_ellipse = dynamic_cast<const Ellipse&>(another);
            if (fabs (another_ellipse.eccentricity() - this->eccentricity()) < EPS
                && fabs(lenght(another_ellipse.f.first, another_ellipse.f.second) - lenght(f.first, f.second)) < EPS) {
                return true;
            } else {
                return false;
            }
        } catch(...) {
            return false;
        }
    }

    virtual bool isSimilarTo(const Shape& another) const override {
        try {
            if (fabs(dynamic_cast<const Ellipse&>(another).eccentricity() - this->eccentricity()) < EPS) {
                return true;
            } else {
                return false;
            }
        } catch(...) {
            return false;
        }
    }

    virtual bool containsPoint(Point point) const override {
        return lenght(f.first, point) + lenght(f.second, point) - 2 * a < EPS;
    }

    virtual void rotate(Point center, double angle) override {
        Point vec (f.first.x - center.x, f.first.y - center.y);
        f.first = Point(center.x + vec.x * cos(angle) - vec.y * sin(angle), center.y + vec.x * sin(angle) + vec.y * cos(angle));

        vec = Point(f.second.x - center.x, f.second.y - center.y);
        f.second = Point(center.x + vec.x * cos(angle) - vec.y * sin(angle), center.y + vec.x * sin(angle) + vec.y * cos(angle));
    }

    virtual void reflex(Point center) override {
        f.first.reflex(center);
        f.second.reflex(center);
    }

    virtual void reflex(Line axis) override {
        f.first = l_reflex(axis, f.first);
        f.second = l_reflex(axis, f.second);
    }

    virtual void scale(Point center, double coefficient) override {
        f.first = Point(center.x + (f.first.x - center.x) * coefficient, center.y + (f.first.y - center.y) * coefficient);
        f.second = Point(center.x + (f.second.x - center.x) * coefficient, center.y + (f.second.y - center.y) * coefficient);
    }

    virtual bool operator==(const Shape& another) const override {
        try {
            Ellipse another_ellipse = dynamic_cast<const Ellipse&>(another);
            if (another_ellipse.focuses() == f && fabs(another_ellipse.a - a) < EPS) {
                return true;
            } else {
                return false;
            }
        } catch(...) {
            return false;
        }
    }

    virtual bool operator!=(const Shape& another) const override {
        return !(*this == another);
    }
};

class Circle: public Ellipse {
protected:
    Point center_point;
    double r = 0.0;

public:
    Circle() = default;

    Circle(const Point center, double r): Ellipse(center, center, 2.0 * r), center_point(center), r(r) {}

    double radius() const {
        return r;
    }

    Point center() const {
        return center_point;
    }

    double perimeter() const override {
        return 2.0 * M_PI * r;
    }

    double area() const override {
        return M_PI * r * r;
    }

    bool containsPoint(Point point) const override {
        return (r - lenght(center_point, point)) > -EPS;
    }

    void rotate(Point center, double angle) override {
        Point vec (this->center_point.x - center.x, this->center_point.y - center.y);
        this->center_point = Point(center.x + vec.x * cos(angle) - vec.y * sin(angle), center.y + vec.x * sin(angle) + vec.y * cos(angle));
    }

    void reflex(Point center) override {
        this->center_point.reflex(center);
    }

    void reflex(Line axis) override {
        this->center_point = l_reflex(axis, center_point);
    }

    void scale(Point center, double coefficient) override {
        r *= coefficient;
        this->center_point = Point(center.x + (this->center_point.x - center.x) * coefficient, center.y + (this->center_point.y - center.y) * coefficient);
    }
};

class Polygon: public Shape {
protected:
    std::vector<Point> vertices;

    //mode = 0 - конгруэнтность, mode = 1 - подобие
    bool isCongruent_isSimilar_Helper(const Polygon& copy, std::vector<Point>& another_vertices, char mode) const {
        Point st = vertices[1] - vertices[0];
        for (size_t i = 0; i < vertices.size(); ++i) {
            Point vec = Point(0.0, 0.0) - another_vertices[i];

            for (size_t j = 0; j < vertices.size(); ++j) {
                another_vertices[j] += vec;
            }

            double angle;
            if (i == vertices.size() - 1) {
                angle = acos(st * (another_vertices[0] - another_vertices[vertices.size() - 1]) / lenght(st, Point(0.0, 0.0))
                            / lenght(another_vertices[0] - another_vertices[vertices.size() - 1], Point(0.0, 0.0)));
            } else {
                angle = acos(st * (another_vertices[i + 1] - another_vertices[i]) / lenght(st, Point(0.0, 0.0))
                            / lenght(another_vertices[i + 1] - another_vertices[i], Point(0.0, 0.0)));
            }

            Polygon t_another (another_vertices);
            t_another.rotate(Point(0.0, 0.0), angle);

            if (mode == 1) {
                double coefficient;
                if (i != vertices.size() - 1) {
                    coefficient = lenght(st, Point(0.0, 0.0)) / lenght(another_vertices[i + 1] - another_vertices[i], Point(0.0, 0.0));
                } else {
                    coefficient = lenght(st, Point(0.0, 0.0)) / lenght(another_vertices[i] - another_vertices[0], Point(0.0, 0.0));
                }
                t_another.scale(Point(0.0, 0.0), coefficient);
            }

            if (copy == t_another) {
                return true;
            }

            t_another.reflex(Line(Point(0.0, 0.0), vertices[i + 1]));
            if (copy == t_another) {
                return true;
            }
            t_another.reflex(Line(Point(0.0, 0.0), vertices[i + 1]));

            t_another.rotate(Point(0.0, 0.0), -angle);
            t_another.rotate(Point(0.0, 0.0), -angle);

            if (copy == t_another) {
                return true;
            }

            for (size_t j = 0; j < vertices.size(); ++j) {
                another_vertices[j] -= vec;
            }
        }
        return false;
    }

public:
    Polygon() = default;

    Polygon(const std::vector<Point>& src_vertices) {
        vertices = src_vertices;
    }

    Polygon(std::initializer_list<Point> lst) {
        vertices = lst;
    }

    virtual size_t verticesCount() const {
        return vertices.size();
    }

    virtual std::vector<Point> getVertices() const {
        return vertices;
    }

    bool isConvex() const {
        bool ans = true;

        int prev_sign = 0;
        for (size_t i = 1; i < vertices.size() - 1; ++i) {
            double product = (vertices[i].x - vertices[i - 1].x) * (vertices[i + 1].y - vertices[i].y)
                           - (vertices[i].y - vertices[i - 1].y) * (vertices[i + 1].x - vertices[i].x);

            if (prev_sign != 0) {
                if ((prev_sign == -1 && product > -EPS) || (prev_sign == 1 && product < -EPS)) {
                    ans = false;
                    break;
                }
            }

            prev_sign = (product > -EPS) ? 1 : -1;
        }

        return ans;
    }

    virtual double perimeter() const override {
        double ans = 0.0;

        for (size_t i = 0; i < vertices.size() - 1; ++i) {
            ans += lenght(vertices[i], vertices[i + 1]);
        }
        ans += lenght(vertices[0], vertices[vertices.size() - 1]);

        return ans;
    }

    virtual double area() const override {
        if (vertices.size() == 0) {
            return 0.0;
        }

        double ans = 0;

        for (size_t i = 0; i < vertices.size() - 1; ++i) {
            ans += vertices[i].y * vertices[i + 1].x - vertices[i].x * vertices[i + 1].y;
        }

        ans += vertices[0].x * vertices[vertices.size() - 1].y - vertices[0].y * vertices[vertices.size() - 1].x;

        return fabs(ans) / 2;
    }

    virtual bool isCongruentTo(const Shape& another) const override {
        try {
            std::vector<Point> another_vertices = dynamic_cast<const Polygon&>(another).getVertices();
            if (another_vertices.size() != vertices.size()) {
                return false;
            } else {
                Point vec = Point(0.0, 0.0) - vertices[0];
                Polygon copy = *this;

                for (size_t i = 0; i < vertices.size(); ++i) {
                    copy.vertices[i] += vec;
                }

                if (isCongruent_isSimilar_Helper(copy, another_vertices, 0)) {
                    return true;
                } else {
                    std::vector<Point> reverse_another_vertices;
                    for (size_t i = 0; i < vertices.size(); ++i) {
                        reverse_another_vertices.push_back(another_vertices[vertices.size() - 1 - i]);
                    }
                    return isCongruent_isSimilar_Helper(copy, reverse_another_vertices, 0);
                }
            }
        return false;
        } catch(...) {
            return false;
        }
    }

    virtual bool isSimilarTo(const Shape& another) const override {
        try {
            std::vector<Point> another_vertices = dynamic_cast<const Polygon&>(another).getVertices();
            if (another_vertices.size() != vertices.size()) {
                return false;
            } else {
                Point vec = Point(0.0, 0.0) - vertices[0];
                Polygon copy = *this;

                for (size_t i = 0; i < vertices.size(); ++i) {
                    copy.vertices[i] += vec;
                }

                if (isCongruent_isSimilar_Helper(copy, another_vertices, 1)) {
                    return true;
                } else {
                    std::vector<Point> reverse_another_vertices;
                    for (size_t i = 0; i < vertices.size(); ++i) {
                        reverse_another_vertices.push_back(another_vertices[vertices.size() - 1 - i]);
                    }
                    return isCongruent_isSimilar_Helper(copy, reverse_another_vertices, 1);
                }
            }
        return false;
        } catch(...) {
            return false;
        }
    }

    virtual bool containsPoint(Point point) const override {
        bool ans = false;

        size_t j = vertices.size() - 1;
        for (size_t i = 0; i < vertices.size(); ++i) {
            if (((vertices[i].y - point.y < -EPS && vertices[j].y - point.y > -EPS) || (vertices[j].y - point.y < -EPS && vertices[i].y - point.y > -EPS)) &&
                (vertices[i].x + (point.y - vertices[i].y) / (vertices[j].y - vertices[i].y) * (vertices[j].x - vertices[i].x) - point.x < -EPS)) {
                ans = !ans;
            }
            j = i;
        }

        return ans;
    }

    virtual void rotate(Point center, double angle) override {
        for (size_t i = 0; i < vertices.size(); ++i) {
            Point vec (vertices[i].x - center.x, vertices[i].y - center.y);
            vertices[i] = Point(center.x + vec.x * cos(angle) - vec.y * sin(angle), center.y + vec.x * sin(angle) + vec.y * cos(angle));
        }
    }

    virtual void reflex(Point center) override {
        for (size_t i = 0; i < vertices.size(); ++i) {
            vertices[i].reflex(center);
        }
    }

    virtual void reflex(Line axis) override {
        for (size_t i = 0; i < vertices.size(); ++i) {
            vertices[i] = l_reflex(axis, vertices[i]);
        }
    }

    virtual void scale(Point center, double coefficient) override {
        for (size_t i = 0; i < vertices.size(); ++i) {
            vertices[i] = Point(center.x + (vertices[i].x - center.x) * coefficient, center.y + (vertices[i].y - center.y) * coefficient);
        }
    }

    virtual bool operator==(const Shape& another) const override {
        try {
            std::vector<Point> another_vertices = dynamic_cast<const Polygon&>(another).getVertices();

            if (another_vertices.size() != vertices.size()) {
                return false;
            }

            size_t start_pos = 0, size = vertices.size();
            while (start_pos < size && vertices[0] != another_vertices[start_pos]) {
                ++start_pos;
            }

            if (start_pos == vertices.size()) {
                return false;
            } else {
                if (another_vertices[(start_pos + 1) % size] != vertices[1] && another_vertices[(start_pos - 1 + size) % size] != vertices[1]) {
                    return false;
                } else {
                    int dir = (another_vertices[(start_pos + 1) % size] == vertices[1]) ? 1 : -1;

                    for (size_t i = (start_pos + 1) % size; i != start_pos; i = (i + dir + size) % size) {
                        if (another_vertices[i] != vertices[(i - start_pos + size) % size]) {
                            return false;
                        }
                    }

                    return true;
                }
            }
        } catch(...) {
            return false;
        }
    }

    virtual bool operator!=(const Shape& another) const override {
        return !(*this == another);
    }
};

class Rectangle: public Polygon {
public:
    Rectangle() = default;

    Rectangle(const Point point_1, const Point point_2, double relation) {
        vertices.push_back(point_1);

        double c = lenght(point_1, point_2);
        double a = c / sqrt(1 + relation * relation);
        double b = sqrt(c * c - a * a);

        double v1_x = (point_2.x - point_1.x) / c;
        double v1_y = (point_2.y - point_1.y) / c;

        if (a < b) {
            std::swap(a, b);
        }

        double v2_x = v1_x * cos(atan(a / b)) - v1_y * sin(atan(a / b));
        double v2_y = v1_x * sin(atan(a / b)) + v1_y * cos(atan(a / b));

        Point t_point (point_1.x + b * v2_x, point_1.y + b * v2_y);
        vertices.push_back(t_point);
        vertices.push_back(point_2);

        std::swap(a, b);

        v2_x = v1_x * cos(atan(a / b)) + v1_y * sin(atan(a / b));
        v2_y = -v1_x * sin(atan(a / b)) + v1_y * cos(atan(a / b));

        t_point = Point(point_1.x + b * v2_x, point_1.y + b * v2_y);
        vertices.push_back(t_point);
    }

    std::pair<Line, Line> diagonals() const {
        return std::make_pair(Line(vertices[0], vertices[2]), Line(vertices[1], vertices[3]));
    }

    Point center() const {
        return middle_point(vertices[0], vertices[2]);
    }

    virtual double perimeter() const override {
        return 2.0 * (lenght(vertices[0], vertices[1]) + lenght(vertices[1], vertices[2]));
    }

    virtual double area() const override {
        if (vertices.size() == 0) {
            return 0.0;
        }

        return lenght(vertices[0], vertices[1]) * lenght(vertices[1], vertices[2]);
    }

    virtual bool containsPoint(Point point) const override {
        double p1 = product(vertices[0], vertices[1], point);
        double p2 = product(vertices[1], vertices[2], point);
        double p3 = product(vertices[2], vertices[3], point);
        double p4 = product(vertices[3], vertices[0], point);

        if ((p1 < -EPS && p2 < -EPS && p3 < -EPS && p4 < -EPS) || (p1 > -EPS && p2 > -EPS && p3 > -EPS && p4 > -EPS)) {
            return true;
        } else {
            return false;
        }
    }

    virtual void rotate(Point center, double angle) override {
        for (size_t i = 0; i < 4; ++i) {
            Point vec (vertices[i].x - center.x, vertices[i].y - center.y);
            vertices[i] = Point(center.x + vec.x * cos(angle) - vec.y * sin(angle), center.y + vec.x * sin(angle) + vec.y * cos(angle));
        }
    }

    virtual void reflex(Point center) override {
        for (size_t i = 0; i < 4; ++i) {
            vertices[i].reflex(center);
        }
    }

    virtual void reflex(Line axis) override {
        for (size_t i = 0; i < 4; ++i) {
            vertices[i] = l_reflex(axis, vertices[i]);
        }
    }

    virtual void scale(Point center, double coefficient) override {
        for (size_t i = 0; i < 4; ++i) {
            vertices[i] = Point(center.x + (vertices[i].x - center.x) * coefficient, center.y + (vertices[i].y - center.y) * coefficient);
        }
    }
};

class Square: public Rectangle {
public:
    Square() = default;

    Square(const Point point_1, const Point point_2) {
        vertices.push_back(point_1);

        double c = lenght(point_1, point_2);
        double a = c / sqrt(2);
        double b = a;

        double v1_x = (point_2.x - point_1.x) / c;
        double v1_y = (point_2.y - point_1.y) / c;

        double v2_x = v1_x * cos(atan(1)) - v1_y * sin(atan(1));
        double v2_y = v1_x * sin(atan(1)) + v1_y * cos(atan(1));

        Point t_point (point_1.x + b * v2_x, point_1.y + b * v2_y);
        vertices.push_back(t_point);
        vertices.push_back(point_2);

        v2_x = v1_x * cos(atan(1)) + v1_y * sin(atan(1));
        v2_y = -v1_x * sin(atan(1)) + v1_y * cos(atan(1));

        t_point = Point(point_1.x + b * v2_x, point_1.y + b * v2_y);
        vertices.push_back(t_point);
    }

    double perimeter() const override {
        return 4.0 * lenght(vertices[0], vertices[1]);
    }

    double area() const override {
        if (vertices.size() == 0) {
            return 0.0;
        }

        return lenght(vertices[0], vertices[1]) * lenght(vertices[0], vertices[1]);
    }

    Circle inscribedCircle() {
        return Circle(this->center(), lenght(vertices[0], vertices[1]) / 2.0);
    }

    Circle circumscribedCircle() {
        return Circle(this->center(), lenght(vertices[0], vertices[2]) / 2.0);
    }
};

class Triangle: public Polygon {
public:
    Triangle() = default;

    Triangle(const Point point_1, const Point point_2, const Point point_3) {// суперконструктор от инициализационного листа?
        vertices.push_back(point_1);
        vertices.push_back(point_2);
        vertices.push_back(point_3);
    }

    double perimeter() const override {
        return lenght(vertices[0], vertices[1]) + lenght(vertices[1], vertices[2]) + lenght(vertices[2], vertices[0]);
    }

    double area() const override {
        if (vertices.size() == 0) {
            return 0.0;
        }

        double a = lenght(vertices[0], vertices[1]), b = lenght(vertices[1], vertices[2]), c = lenght(vertices[0], vertices[2]);
        double p = (a + b + c) / 2.0;

        return sqrt(p * (p - a) * (p - b) * (p - c));
    }

    Point orthocenter() const {
        Point vec1 = get_vec_from_2_points(vertices[0], vertices[1]);
        vec1 = Point(-vec1.y, vec1.x);

        Point vec2 = get_vec_from_2_points(vertices[1], vertices[2]);;
        vec2 = Point(-vec2.y, vec2.x);

        Line h1 (vertices[2], vertices[2] + vec1);
        Line h2 (vertices[0], vertices[0] + vec2);

        double det = h1.a * h2.b - h2.a * h1.b;

        return Point((h1.a * (-h2.c) - h2.a * (-h1.c)) / det, ((-h1.c) * h2.b - (-h2.c) * h1.b) / det);
    }

    Point centroid() const {
        return Point((vertices[0].x + vertices[1].x + vertices[2].x) / 3, (vertices[0].y + vertices[1].y + vertices[2].y) / 3);
    }

    Point circumscribedCircleCenter() const {
        Point vec1 = get_vec_from_2_points(vertices[0], vertices[1]);
        vec1 = Point(-vec1.y, vec1.x);

        Point vec2 = get_vec_from_2_points(vertices[1], vertices[2]);
        vec2 = Point(-vec2.y, vec2.x);

        Line mp1 (middle_point(vertices[0], vertices[1]), middle_point(vertices[0], vertices[1]) + vec1);
        Line mp2 (middle_point(vertices[1], vertices[2]), middle_point(vertices[1], vertices[2]) + vec2);

        double det = mp1.a * mp2.b - mp2.a * mp1.b;

        return Point((mp1.a * (-mp2.c) - mp2.a * (-mp1.c)) / det, ((-mp1.c) * mp2.b - (-mp2.c) * mp1.b) / det);
    }

    Point inscribedCircleCenter() const {
        Point vec1 = get_vec_from_2_points(vertices[0], vertices[1]);
        Point vec2 = get_vec_from_2_points(vertices[1], vertices[0]);
        Point vec3 = get_vec_from_2_points(vertices[1], vertices[2]);
        Point vec4 = get_vec_from_2_points(vertices[0], vertices[2]);

        Point t_point_1 = vertices[0] + vec1;
        Point t_point_2 = vertices[0] + vec4;
        Line b1 (vertices[0], middle_point(t_point_1, t_point_2));

        t_point_1 = vertices[1] + vec2;
        t_point_2 = vertices[1] + vec3;
        Line b2 (vertices[1], middle_point(t_point_1, t_point_2));

        double det = b1.a * b2.b - b2.a * b1.b;

        return Point((b1.a * (-b2.c) - b2.a * (-b1.c)) / det, ((-b1.c) * b2.b - (-b2.c) * b1.b) / det);
    }

    Circle circumscribedCircle() const {
        return Circle(this->circumscribedCircleCenter(),
                      lenght(vertices[0], vertices[1]) * lenght(vertices[1], vertices[2]) * lenght(vertices[0], vertices[2]) / (4.0 * this->area()));
    }

    Circle inscribedCircle() const {
        return Circle(this->inscribedCircleCenter(), 2.0 * this->area() / this->perimeter());
    }

    Line EulerLine() {
        if (this->centroid() == this->orthocenter()) {
            return Line(this->centroid(), 0.0);
        } else {
            return Line(this->orthocenter(), this->centroid());
        }
    }

    Circle ninePointsCircle() {
        Point oc = this->orthocenter();
        Point cc = this->circumscribedCircleCenter();
        return Circle(Point((oc.x + cc.x) / 2.0, (oc.y + cc.y) / 2.0),
                      lenght(vertices[0], vertices[1]) * lenght(vertices[1], vertices[2]) * lenght(vertices[0], vertices[2]) / (8.0 * this->area()));
    }

    bool containsPoint(Point point) const override {
        double p1 = product(vertices[0], vertices[1], point);
        double p2 = product(vertices[1], vertices[2], point);
        double p3 = product(vertices[2], vertices[0], point);

        if ((p1 < -EPS && p2 < -EPS && p3 < -EPS) || (p1 > -EPS && p2 > -EPS && p3 > -EPS)) {
            return true;
        } else {
            return false;
        }
    }

    void rotate(Point center, double angle) override {
        for (size_t i = 0; i < 3; ++i) {
            Point vec (vertices[i].x - center.x, vertices[i].y - center.y);
            vertices[i] = Point(center.x + vec.x * cos(angle) - vec.y * sin(angle), center.y + vec.x * sin(angle) + vec.y * cos(angle));
        }
    }

    void reflex(Point center) override {
        for (size_t i = 0; i < 3; ++i) {
            vertices[i].reflex(center);
        }
    }

    void reflex(Line axis) override {
        for (size_t i = 0; i < 3; ++i) {
            vertices[i] = l_reflex(axis, vertices[i]);
        }
    }

    void scale(Point center, double coefficient) override {
        for (size_t i = 0; i < 3; ++i) {
            vertices[i] = Point(center.x + (vertices[i].x - center.x) * coefficient, center.y + (vertices[i].y - center.y) * coefficient);
        }
    }
};
