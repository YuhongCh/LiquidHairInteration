#pragma once

struct Color {
	static constexpr Scalar inv255 = 1.0 / 255.0;

	union {
		struct {
			Integer r;
			Integer g;
			Integer b;
			Integer a;
		};
		Integer rgba[4];
	};

	Color(const Integer& _r = 0, const Integer& _g = 0, const Integer& _b = 0, const Integer& _a = 0)
		: r(_r), g(_g), b(_b), a(_a) {
		if (r < 0 || r > 255 || g < 0 || g > 255 || b < 0 || b > 255 || a < 0 || a > 255) {
			throw std::invalid_argument("rgba value must be in range [0, 255]");
			return;
		}
	}

	inline std::tuple<Scalar, Scalar, Scalar, Scalar> toScalar() const {
		return { r * inv255, g * inv255, b * inv255, a * inv255 };
	}

	inline static Color Red() { return Color(255, 0, 0, 255); }
	inline static Color Green() { return Color(0, 255, 0, 255); }
	inline static Color Blue() { return Color(0, 0, 255, 255); }
	inline static Color White() { return Color(255, 255, 255, 255); }
	inline static Color Black() { return Color(0, 0, 0, 255); }
	inline static Color Gray() { return Color(128, 128, 128, 255); }
};
