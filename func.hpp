#ifndef FUNC_H

#define FUNC_H

#include <functional>

template <typename out, typename in>
out e(in value = in{}) {
    return out{};
}

template <typename out = double, typename in = double>
class func {
public:

    func() : func(e<out, in>){}

    func (const std::function<out(in)>& f_) : f(f_) {}

    func(out (*p)(in)) : func(std::function<out(in)>(p)) {}

    func(const func<in, out>& f_) : func(f_.f ) {}

    func(const in& value) {
        *this = constant(value);
    }

    out operator()(in value = in{}) const {
        return f(value);
    }

    // Stores given constant in f
    func& constant(const in& value) {
        this->f = [value](in vd) -> out { return value; };
        return *this;
    }

    // Returns given value (may convert from in_type to out_type)
    func& chi() {
        this-> f = [](in value) -> out { return value; };
        return *this;
    }

    func& operator=(const func& f_) {
        this->f = f_.f;
        return *this;
    }

    func& operator+=(const func& f_) {
        this->f = [*this, f_](in value) -> out { return f(value) + f_(value); };
        return *this;
    }

    func& operator-=(const func& f_) {
        this->f = [*this, f_](in value) -> out { return f(value) - f_(value); };
        return *this;
    }

    func operator-() const {
        return func() - *this;
    }

    func& operator*=(const func& f_) {
        this->f = [*this, f_](in value) -> out { return f(value) * f_(value); };
        return *this;
    }

    func& operator/=(const func& f_) {
        this->f = [*this, f_](in value) -> out { return f(value) / f_(value); };
        return *this;
    }

    func& pow_to(const func& f_) {
        this->f = [*this, f_](in value) -> out { return std::pow(f(value), f_(value)); };
        return *this;
    }

    func pow(const func& f_) const {
        func<out, in> f__ = *this;
        return pow_to(f__, f_);
    }

    func& compose_to(const func& f_) {
        this->f = [*this, f_](in value) -> out { return f(f_(value)); };
        return *this;
    }

    func compose(const func& f_) {
        func<out, in> f__ = *this;
        return compose_to(f__, f_);        
    }

    bool operator==(const func& f_) {
        return false;
    }

private:
    std::function<out(in)> f;

};

template <class out, class in>
std::ostream& operator<<(std::ostream& ost, const func<out, in>& f) {
    ost << f(in{});
    return ost;
}

template <class out, class in>
std::istream& operator>>(std::istream& ist, func<out, in>& f) {
    in value;
    ist >> value;
    f.constant(value);
    return ist;
}

template<class out, class in>
inline func<out, in> operator+(const func<out, in>& f1, const func<out, in>& f2) {
    func<out, in> f3 = f1;
    return f3 += f2;
}

template<class out, class in>
inline func<out, in> operator-(const func<out, in>& f1, const func<out, in>& f2) {
    func<out, in> f3 = f1;
    return f3 -= f2;
}

template<class out, class in>
inline func<out, in> operator*(const func<out, in>& f1, const func<out, in>& f2) {
    func<out, in> f3 = f1;
    return f3 *= f2;
}

template<class out, class in>
inline func<out, in> operator/(const func<out, in>& f1, const func<out, in>& f2) {
    func<out, in> f3 = f1;
    return f3 /= f2;
}


#endif