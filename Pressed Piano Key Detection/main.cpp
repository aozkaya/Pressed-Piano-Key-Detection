#include <opencv2\imgproc\imgproc.hpp>
#include <opencv2\opencv.hpp>
#include <opencv2\highgui\highgui.hpp>
#include <opencv2\objdetect\objdetect.hpp>
#include <stdio.h>
#include <iostream>
#include "midi.h"
#include <Windows.h>
#include <MMSystem.h>
//top filter
int Y_MIN = 70;
int Y_MAX = 198;
int Cr_MIN = 141;
int Cr_MAX = 255;
int Cb_MIN = 0;
int Cb_MAX = 256;
//front BGR filter
int B_MIN = 1;
int B_MAX = 255;
int G_MIN = 1;
int G_MAX = 255;
int R_MIN = 1;
int R_MAX = 255;
//front YCrCb filter
int Y2_MIN = 92;
int Y2_MAX = 217;
int Cr2_MIN = 0;
int Cr2_MAX = 180;
int Cb2_MIN = 0;
int Cb2_MAX = 255;
//front HLS filter
int H_MIN = 0;
int H_MAX = 105;
int L_MIN = 0;
int L_MAX = 255;
int S_MIN = 0;
int S_MAX = 255;
//top YCrCb key filter
int K_Y_MIN = 100;
int K_Y_MAX = 255;
int K_Cr_MIN = 0;
int K_Cr_MAX = 145;
int K_Cb_MIN = 0;
int K_Cb_MAX = 255;
//top YCrCb red filter
int R_Y_MIN = 0;
int R_Y_MAX = 255;
int R_Cr_MIN = 129;
int R_Cr_MAX = 159;
int R_Cb_MIN = 0;
int R_Cb_MAX = 255;
//top HLS red filter
int R_H_MIN = 0;
int R_H_MAX = 255;
int R_L_MIN = 0;
int R_L_MAX = 95;
int R_S_MIN = 50;
int R_S_MAX = 255;
//front YCrCb key filter
int K_Y2_MIN = 0;
int K_Y2_MAX = 255;
int K_Cr2_MIN = 0;
int K_Cr2_MAX = 255;
int K_Cb2_MIN = 117;
int K_Cb2_MAX = 130;
//front HLS key filter
int K_H_MIN = 0;
int K_H_MAX = 255;
int K_L_MIN = 145;
int K_L_MAX = 255;
int K_S_MIN = 0;
int K_S_MAX = 255;
//world info
int camera_width = 1920;
int camera_height = 1080;
float distUpper = 43.5;
float realKeyLength = 14.6;
float zUpper = 4350;

using namespace cv;
using namespace std;




void on_trackbar(int, void*)
{

}


float distance(Point x, Point y) {
	return sqrt((x.x - y.x) * (x.x - y.x) + (x.y - y.y) * (x.y - y.y));
}

float angle(Point x, Point y, Point z) {
	float dot = (x.x - y.x) * (z.x - y.x) + (x.y - y.y) * (z.y - y.y);
	float angle = acos(dot / (distance(x, y) * distance(y, z) * 2));
	return angle * 180 / CV_PI;
}

Point getIntersection(Point p1, Point p2, Point x, int mode) {
	if (mode == 1) {
		if (p1.x - p2.x == 0) return x;
		else {
			int y = p2.y * abs(p1.x - x.x) / abs(p1.x - p2.x) + p1.y * abs(p2.x - x.x) / abs(p1.x - p2.x);
			return Point(x.x, y);
		}
	}
	else {
		if (p1.y - p2.y == 0) return x;
		else {
			int xx = p2.x * abs(p1.y - x.y) / abs(p1.y - p2.y) + p1.x * abs(p2.y - x.y) / abs(p1.y - p2.y);
			return Point(xx, x.y);
		}
	}
}

float slope(Point pt1, Point pt2) {
	if (pt2.x - pt1.x == 0) return 100000000000;
	else return (1.0 * pt2.y - pt1.y) / (1.0 * pt2.x - pt1.x);
}


bool checkoverlap(Point x, vector<Point> y, int num) {
	for (int i = 0; i < num; i++) {
		if (distance(y[i], x) < 40) return false;
	}
	return true;
}

float distancePointLine(Point p, Vec4f l) {
	Point pointA = Point(l[0], l[1]);
	Point pointB = Point(l[2], l[3]);

	int A = 0, B = 0, C = 0;
	A = pointA.y - pointB.y;
	B = pointB.x - pointA.x;
	C = pointA.x * pointB.y - pointA.y * pointB.x;

	float distance = 0;
	distance = ((float)abs(A * p.x + B * p.y + C)) / ((float)sqrtf(A * A + B * B));
	return abs(distance);
}

bool sortHelper(const Point& a, const Point& b) {
	if (a.x < b.x) return true;
	else return false;
}

bool sortHelper2(const Vec4i& a, const Vec4i& b) {
	if (a[0] < b[0]) return true;
	else return false;
}

bool sortHelper3(const Point& a, const Point& b) {
	if (a.y > b.y) return true;
	else return false;
}

bool sortHelper4(const int& a, const int& b) {
	if (a < b) return true;
	else return false;
}


void createTrackbars() {
	namedWindow("Top Trackbars", 0);

	namedWindow("Front Trackbars", 0);

	createTrackbar("Y_MIN", "Top Trackbars", &Y_MIN, 255, on_trackbar);
	createTrackbar("Y_MAX", "Top Trackbars", &Y_MAX, 255, on_trackbar);
	createTrackbar("Cr_MIN", "Top Trackbars", &Cr_MIN, 255, on_trackbar);
	createTrackbar("Cr_MAX", "Top Trackbars", &Cr_MAX, 255, on_trackbar);
	createTrackbar("Cb_MIN", "Top Trackbars", &Cb_MIN, 255, on_trackbar);
	createTrackbar("Cb_MAX", "Top Trackbars", &Cb_MAX, 255, on_trackbar);

	createTrackbar("B_MIN", "Front Trackbars", &B_MIN, 255, on_trackbar);
	createTrackbar("B_MAX", "Front Trackbars", &B_MAX, 255, on_trackbar);
	createTrackbar("G_MIN", "Front Trackbars", &G_MIN, 255, on_trackbar);
	createTrackbar("G_MAX", "Front Trackbars", &G_MAX, 255, on_trackbar);
	createTrackbar("R_MIN", "Front Trackbars", &R_MIN, 255, on_trackbar);
	createTrackbar("R_MAX", "Front Trackbars", &R_MAX, 255, on_trackbar);

	createTrackbar("Y_MIN", "Front Trackbars", &Y2_MIN, 255, on_trackbar);
	createTrackbar("Y_MAX", "Front Trackbars", &Y2_MAX, 255, on_trackbar);
	createTrackbar("Cr_MIN", "Front Trackbars", &Cr2_MIN, 255, on_trackbar);
	createTrackbar("Cr_MAX", "Front Trackbars", &Cr2_MAX, 255, on_trackbar);
	createTrackbar("Cb_MIN", "Front Trackbars", &Cb2_MIN, 255, on_trackbar);
	createTrackbar("Cb_MAX", "Front Trackbars", &Cb2_MAX, 255, on_trackbar);

	createTrackbar("H_MIN", "Front Trackbars", &H_MIN, 255, on_trackbar);
	createTrackbar("H_MAX", "Front Trackbars", &H_MAX, 255, on_trackbar);
	createTrackbar("L_MIN", "Front Trackbars", &L_MIN, 255, on_trackbar);
	createTrackbar("L_MAX", "Front Trackbars", &L_MAX, 255, on_trackbar);
	createTrackbar("S_MIN", "Front Trackbars", &S_MIN, 255, on_trackbar);
	createTrackbar("S_MAX", "Front Trackbars", &S_MAX, 255, on_trackbar);
}

void calculateFingers(Mat canny, Mat imageCut, vector<Point>& fingertips, int mode, vector<Point> whiteLine) {
	vector<vector<Point>> contours;
	vector<Point> Rfingertips;
	int handContour = 0;
	int hand2Contour = 0;
	findContours(canny, contours, RETR_TREE, CHAIN_APPROX_SIMPLE);
	for (int i = 0; i < contours.size(); i++) {
		if (contours[i].size() > contours[handContour].size()) { hand2Contour = handContour;  handContour = i; }
		else if (contours[i].size() > contours[hand2Contour].size() || hand2Contour == handContour) { hand2Contour = i; }
	}

	if (contours.size() < 2) { cout << "at least one of the hands cannot be detected" << endl; }
	else {


		vector<vector<int>> indices = vector<vector<int>>(contours.size());
		vector<vector<Vec4i>> defects = vector<vector<Vec4i>>(contours.size());

		convexHull(contours[handContour], indices[handContour], true, false);
		convexHull(contours[hand2Contour], indices[hand2Contour], true, false);

		sort(indices[handContour].begin(), indices[handContour].end(), sortHelper4);
		sort(indices[hand2Contour].begin(), indices[hand2Contour].end(), sortHelper4);

		fingertips = vector<Point>(100);
		Rfingertips = vector<Point>(100);
		int ctL = 0;
		int ctR = 0;

		if (indices[handContour].size() > 3) {
			convexityDefects(contours[handContour], indices[handContour], defects[handContour]);
			for (int i = 0; i < defects[handContour].size(); i++) {
				Vec4i vec = defects[handContour][i];
				Point Point1(contours[handContour][vec[0]]);
				Point Point2(contours[handContour][vec[1]]);
				Point Point3(contours[handContour][vec[2]]);

				if (angle(Point1, Point3, Point2) < 90)
				{
					if (checkoverlap(Point2, Rfingertips, ctR) && mode * (whiteLine[0].y - Point2.y) > 0) {
						fingertips[ctL] = Point2;
						ctL++;
					}
					if (checkoverlap(Point1, fingertips, ctL) && mode * (whiteLine[0].y - Point1.y) > 0) {
						Rfingertips[ctR] = Point1;
						ctR++;
					}
				}
			}
		}
		if (indices[hand2Contour].size() > 3) {
			convexityDefects(contours[hand2Contour], indices[hand2Contour], defects[hand2Contour]);
			for (int i = 0; i < defects[hand2Contour].size(); i++) {
				Vec4i vec = defects[hand2Contour][i];
				Point Point1(contours[hand2Contour][vec[0]]);
				Point Point2(contours[hand2Contour][vec[1]]);
				Point Point3(contours[hand2Contour][vec[2]]);

				if (angle(Point1, Point3, Point2) < 90)
				{
					if (checkoverlap(Point2, Rfingertips, ctR) && mode * (whiteLine[0].y - Point2.y) > 0) {
						fingertips[ctL] = Point2;
						ctL++;
					}
					if (checkoverlap(Point1, fingertips, ctL) && mode * (whiteLine[0].y - Point1.y) > 0) {
						Rfingertips[ctR] = Point1;
						ctR++;
					}
				}

			}
		}
		for (int i = 0; i < ctR; i++) {
			fingertips[ctL + i] = Rfingertips[i];
		}
		fingertips.resize(ctL + ctR);
	}
}

void calculateFingers2(Mat canny, Mat imageCut, vector<Point>& fingertips, int mode, vector<Point> whiteLine) {
	vector<vector<Point>> contours;
	vector<Point> Rfingertips;
	int convexHullDist;
	if (mode == 1) convexHullDist = 70;
	else convexHullDist = 30;

	int handContour = 0;
	int hand2Contour = 0;
	findContours(canny, contours, RETR_TREE, CHAIN_APPROX_SIMPLE);
	if (contours.size() < 2) {
		cout << "at least one of the hands cannot be detected" << endl;
	}
	else {
		for (int i = 0; i < contours.size(); i++) {
			if (contours[i].size() > contours[handContour].size()) { hand2Contour = handContour;  handContour = i; }
			else if (contours[i].size() > contours[hand2Contour].size() || hand2Contour == handContour) { hand2Contour = i; }
		}

		vector<vector<Point>> points = vector<vector<Point>>(contours.size());

		convexHull(contours[handContour], points[handContour], false, true);
		convexHull(contours[hand2Contour], points[hand2Contour], false, true);

		fingertips = vector<Point>(100);
		int ct = 0;

		int status = 0;
		for (int i = 0; i < contours[handContour].size(); i++) {
			if (status == 0) {
				if (contours[handContour][i].y < contours[handContour][contours[handContour].size() - 1].y) status = -1;
				if (contours[handContour][i].y > contours[handContour][contours[handContour].size() - 1].y) status = 1;
			}
			else if (status == 1) {
				if (contours[handContour][i].y < contours[handContour][i - 1].y) {
					status = -1;
					if (mode == -1) {
						for (int j = 0; j < points[handContour].size() - 1; j++) {
							if ((contours[handContour][i - 1].x <= points[handContour][j].x && contours[handContour][i - 1].x >= points[handContour][j + 1].x) ||
								(contours[handContour][i - 1].x >= points[handContour][j].x && contours[handContour][i - 1].x <= points[handContour][j + 1].x)) {
								Point intersection = getIntersection(points[handContour][j], points[handContour][j + 1], contours[handContour][i - 1], 1);
								if (abs(contours[handContour][i - 1].y - intersection.y) < convexHullDist && mode * (whiteLine[0].y - contours[handContour][i - 1].y) > 0) {
									fingertips[ct] = contours[handContour][i - 1];
									ct++;
									break;
								}
							}
						}
					}
				}
			}
			else {
				if (contours[handContour][i].y > contours[handContour][i - 1].y) {
					status = 1;
					if (mode == 1) {
						for (int j = 0; j < points[handContour].size() - 1; j++) {
							if ((contours[handContour][i - 1].x <= points[handContour][j].x && contours[handContour][i - 1].x >= points[handContour][j + 1].x) ||
								(contours[handContour][i - 1].x >= points[handContour][j].x && contours[handContour][i - 1].x <= points[handContour][j + 1].x)) {
								Point intersection = getIntersection(points[handContour][j], points[handContour][j + 1], contours[handContour][i - 1], 1);
								if (abs(contours[handContour][i - 1].y - intersection.y) < convexHullDist && mode * (whiteLine[0].y - contours[handContour][i - 1].y) > 0) {
									fingertips[ct] = contours[handContour][i - 1];
									ct++;
									break;
								}
							}
						}
					}
				}
			}
		}

		status = 0;
		for (int i = 0; i < contours[hand2Contour].size(); i++) {
			if (status == 0) {
				if (contours[hand2Contour][i].y < contours[hand2Contour][contours[hand2Contour].size() - 1].y) status = -1;
				if (contours[hand2Contour][i].y > contours[hand2Contour][contours[hand2Contour].size() - 1].y) status = 1;
			}
			else if (status == 1) {
				if (contours[hand2Contour][i].y < contours[hand2Contour][i - 1].y) {
					status = -1;
					if (mode == -1) {
						for (int j = 0; j < points[hand2Contour].size() - 1; j++) {
							if ((contours[hand2Contour][i - 1].x <= points[hand2Contour][j].x && contours[hand2Contour][i - 1].x >= points[hand2Contour][j + 1].x) ||
								(contours[hand2Contour][i - 1].x >= points[hand2Contour][j].x && contours[hand2Contour][i - 1].x <= points[hand2Contour][j + 1].x)) {
								Point intersection = getIntersection(points[hand2Contour][j], points[hand2Contour][j + 1], contours[hand2Contour][i - 1], 1);
								if (abs(contours[hand2Contour][i - 1].y - intersection.y) < convexHullDist && mode * (whiteLine[0].y - contours[hand2Contour][i - 1].y) > 0) {
									fingertips[ct] = contours[hand2Contour][i - 1];
									ct++;
									break;
								}
							}
						}
					}
				}
			}
			else {
				if (contours[hand2Contour][i].y > contours[hand2Contour][i - 1].y) {
					status = 1;
					if (mode == 1) {
						for (int j = 0; j < points[hand2Contour].size() - 1; j++) {
							if ((contours[hand2Contour][i - 1].x <= points[hand2Contour][j].x && contours[hand2Contour][i - 1].x >= points[hand2Contour][j + 1].x) ||
								(contours[hand2Contour][i - 1].x >= points[hand2Contour][j].x && contours[hand2Contour][i - 1].x <= points[hand2Contour][j + 1].x)) {
								Point intersection = getIntersection(points[hand2Contour][j], points[hand2Contour][j + 1], contours[hand2Contour][i - 1], 1);
								if (abs(contours[hand2Contour][i - 1].y - intersection.y) < convexHullDist && mode * (whiteLine[0].y - contours[hand2Contour][i - 1].y) > 0) {
									fingertips[ct] = contours[hand2Contour][i - 1];
									ct++;
									break;
								}
							}
						}
					}
				}
			}
		}
		fingertips.resize(ct);
	}
}

void fillMissingFingers(Mat image, Mat imageCut, Mat imageCut2, vector<Point>& fingertips, vector<Point>& fingertips2, vector<Point> mFingertips, vector<Point> mFingertips2,
	vector<int>& topIndex, vector<int>& frontIndex, vector<int> topIndex2, vector<int> frontIndex2, vector<float>& topKeyProportion, vector<float>& frontKeyProportion, vector<float> topKeyProportion2, vector<float> frontKeyProportion2, Mat image2) {

	int hullFingerCount = fingertips.size();
	fingertips.resize(hullFingerCount + mFingertips.size());
	topIndex.resize(hullFingerCount + topIndex2.size());
	topKeyProportion.resize(hullFingerCount + topKeyProportion2.size());
	int ct = 0;

	for (int i = 0; i < mFingertips.size(); i++) {
		int match = 0;
		for (int j = 0; j < hullFingerCount + ct; j++) {
			if (topIndex2[i] == topIndex[j]) {
				match = 1;
				if (fingertips[j].y > mFingertips[i].y) {
					fingertips[j] = mFingertips[i];
					topKeyProportion[j] = topKeyProportion2[i];
				}
			}
			if (distance(mFingertips[i], fingertips[j]) < 30) {
				match = 1;
				if (fingertips[j].y > mFingertips[i].y) {
					fingertips[j] = mFingertips[i];
					topIndex[j] = topIndex2[i];
					topKeyProportion[j] = topKeyProportion2[i];
				}
			}
		}
		if (match == 0) {
			fingertips[hullFingerCount + ct] = mFingertips[i];
			topIndex[hullFingerCount + ct] = topIndex2[i];
			topKeyProportion[hullFingerCount + ct] = topKeyProportion2[i];
			ct++;
		}
	}

	fingertips.resize(hullFingerCount + ct);
	topIndex.resize(hullFingerCount + ct);
	topKeyProportion.resize(hullFingerCount + ct);

	for (int j = 0; j < fingertips.size(); j++) {
		circle(image, fingertips[j], 3, Scalar(255, 100, 255), 3);
	}

	hullFingerCount = fingertips2.size();
	fingertips2.resize(hullFingerCount + mFingertips2.size());
	frontIndex.resize(hullFingerCount + frontIndex2.size());
	frontKeyProportion.resize(hullFingerCount + frontKeyProportion2.size());
	ct = 0;

	for (int i = 0; i < mFingertips2.size(); i++) {
		int match = 0;
		for (int j = 0; j < hullFingerCount + ct; j++) {
			if (distance(mFingertips2[i], fingertips2[j]) < 50) {
				match = 1;
				if (fingertips2[j].y < mFingertips2[i].y) {
					fingertips2[j] = mFingertips2[i];
					frontIndex[j] = frontIndex2[i];
					frontKeyProportion[j] = frontKeyProportion2[i];
				}
			}
		}
		if (match == 0) {
			fingertips2[hullFingerCount + ct] = mFingertips2[i];
			frontIndex[hullFingerCount + ct] = frontIndex2[i];
			frontKeyProportion[hullFingerCount + ct] = frontKeyProportion2[i];
			ct++;
		}
	}

	fingertips2.resize(hullFingerCount + ct);
	frontIndex.resize(hullFingerCount + ct);
	frontKeyProportion.resize(hullFingerCount + ct);

	for (int j = 0; j < fingertips2.size(); j++) {
		circle(image2, fingertips2[j], 3, Scalar(0, 255, 0), 3);
	}
}

Mat calculateCanny(Mat imageCut, Mat imageFiltered, bool display) {
	Mat canny;
	Mat imageBW;
	Mat imageEdges;
	cvtColor(imageCut, imageBW, COLOR_BGR2GRAY);
	Canny(imageBW, imageEdges, 60, 120, 3);
	dilate(imageEdges, imageEdges, getStructuringElement(MORPH_RECT, Size(7, 7), Point(-1, -1)));
	canny = imageFiltered - imageEdges;
	dilate(canny, canny, getStructuringElement(MORPH_RECT, Size(3, 3), Point(-1, -1)));
	return canny;
}

void longLine(Mat cdst, Point pt1, Point pt2, Scalar color) {
	Point start, finish;
	float sl = slope(pt1, pt2);
	if (abs(sl) > 1) {
		start.x = 0;
		finish.x = camera_width;
		start.y = -pt1.x * sl + pt1.y;
		finish.y = -(pt2.x - camera_width) * sl + pt2.y;
		line(cdst, start, finish, color, 2, LINE_AA);
	}
}


Point findIntersection(Vec4i LineA, Vec4i LineB) {
	double ka, kb;
	ka = (1.0 * LineA[3] - LineA[1]) / (1.0 * LineA[2] - LineA[0]); //Find the slope of LineA
	kb = (1.0 * LineB[3] - LineB[1]) / (1.0 * LineB[2] - LineB[0]); //Find the slope of LineB

	Point2f crossPoint;
	crossPoint.x = (ka * LineA[0] - LineA[1] - kb * LineB[0] + LineB[1]) / (ka - kb);
	crossPoint.y = (ka * kb * (1.0 * LineA[0] - LineB[0]) + ka * LineB[1] - kb * LineA[1]) / (ka - kb);
	return crossPoint;
}

void findVanishingPoint(vector<Vec4i> lines, Mat cdst, Point& vanishingPoint) {

	vector<Point> candidatePoints = vector<Point>(lines.size() * lines.size());
	int counter = 0;
	Point intersection;

	for (size_t i = 0; i < lines.size(); i++) {
		for (size_t j = i + 1; j < lines.size(); j++) {

			Vec4i l = lines[i];
			Vec4i l1 = lines[j];

			float sl = slope(Point(l[0], l[1]), Point(l[2], l[3]));
			float sl1 = slope(Point(l1[0], l1[1]), Point(l1[2], l1[3]));
			if (abs(sl) > 1 && abs(sl1) > 1 && (sl1 > 0 && sl < 0 || sl1 < 0 && sl >0)) {

				intersection = findIntersection(l, l1);
				if (intersection.y > 0) {
					candidatePoints[counter] = intersection;
					counter++;
				}
			}
		}
	}
	float xSum = 0;
	float ySum = 0;

	for (size_t i = 0; i < counter; i++) {

		xSum = xSum + candidatePoints[i].x;
		ySum = ySum + candidatePoints[i].y;

	}

	Point avgIntersection = Point(xSum / counter, ySum / counter);
	vanishingPoint = Point(avgIntersection.x, avgIntersection.y - camera_height);
}


void Hough(Mat dst, Point& vanishingPoint, vector<Point>& whiteLine, vector<Point>& blackLine, bool mode) {

	Mat cdst;
	Mat big(dst.rows * 2, dst.cols, dst.type(), Scalar(0, 0, 0));
	dst.copyTo(big(Rect(0, dst.rows, dst.cols, dst.rows)));
	cvtColor(big, cdst, COLOR_GRAY2BGR);

	vector<Vec4i> linesP;
	HoughLinesP(big, linesP, 1, CV_PI / 180, 50, 100, 10); // probabilistic hough line transform
	if (mode) {
		for (size_t i = 0; i < linesP.size(); i++) {
			Vec4i l = linesP[i];
			longLine(cdst, Point(l[0], l[1]), Point(l[2], l[3]), Scalar(150, 0, 0));
		}
	}
	vector<Vec2f> lines;
	if (mode) HoughLines(big, lines, 0.5, CV_PI / 180, 80, 0, 0); // hough line transform
	else HoughLines(big, lines, 0.5, CV_PI / 180, 110, 0, 0);
	vector<int> ylist = vector<int>(lines.size());
	int lasty = 0;
	vector<int> ilist = vector<int>(lines.size());
	int lasti = 0;
	for (size_t i = 0; i < lines.size(); i++)
	{

		float rho = lines[i][0], theta = lines[i][1];
		Point pt1, pt2;
		double a = cos(theta), b = sin(theta);
		double x0 = a * rho, y0 = b * rho;
		pt1.x = cvRound(x0 + camera_width * (-b));
		pt1.y = cvRound(y0 + camera_height * (a));
		pt2.x = cvRound(x0 - camera_width * (-b));
		pt2.y = cvRound(y0 - camera_height * (a));
		bool isvalid = true;
		for (size_t j = 0; j < lasty; j++) {
			if (abs(pt1.y - ylist[j]) < camera_height / 80 || abs(slope(pt1, pt2)) > 0.1) isvalid = false;
		}
		if (isvalid) {
			ylist[lasty] = pt1.y;
			lasty++;
			ilist[lasti] = i;
			lasti++;
		}
	}
	for (size_t i = 0; i < lasti; i++) {
		float rho = lines[ilist[i]][0], theta = lines[ilist[i]][1];
		Point pt1, pt2;
		double a = cos(theta), b = sin(theta);
		double x0 = a * rho, y0 = b * rho;
		pt1.x = cvRound(x0 + camera_width * (-b));
		pt1.y = cvRound(y0 + camera_height * (a));
		pt2.x = cvRound(x0 - camera_width * (-b));
		pt2.y = cvRound(y0 - camera_height * (a));
		bool iswhite = true;
		bool iswhiteT = true;
		int isblack = 0;
		line(cdst, pt1, pt2, Scalar(100, 100, 100), 2, LINE_AA);
		for (size_t j = 0; j < lasty; j++) {
			if (pt1.y > ylist[j]) iswhite = false;
			if (pt1.y < ylist[j]) iswhiteT = false;
			if (pt1.y < ylist[j]) isblack++;
		}
		if (iswhite && mode) {
			line(cdst, pt1, pt2, Scalar(0, 0, 255), 2, LINE_AA);
			pt1.y = pt1.y - camera_height;
			pt2.y = pt2.y - camera_height;
			whiteLine[0] = pt1;
			whiteLine[1] = pt2;
		}
		if (iswhiteT && !mode) {
			line(cdst, pt1, pt2, Scalar(0, 0, 255), 2, LINE_AA);
			pt1.y = pt1.y - camera_height;
			pt2.y = pt2.y - camera_height;
			whiteLine[0] = pt1;
			whiteLine[1] = pt2;
		}
		if (isblack == 1) {
			line(cdst, pt1, pt2, Scalar(0, 255, 0), 2, LINE_AA);
			pt1.y = pt1.y - camera_height;
			pt2.y = pt2.y - camera_height;
			blackLine[0] = pt1;
			blackLine[1] = pt2;
		}
	}
	if (mode) findVanishingPoint(linesP, cdst, vanishingPoint);
}

void upperToWorld(vector<Point> fingertips, Mat imageCut, vector<Point> whiteLine, vector<Point>& worldFingers, int whiteKeyLength) {
	worldFingers = vector<Point>(fingertips.size());
	vector<float> cross = vector<float>(fingertips.size());

	float fUpper = whiteKeyLength * distUpper / realKeyLength;
	for (int i = 0; i < fingertips.size(); i++) {
		worldFingers[i].x = (fingertips[i].x - camera_width / 2) * zUpper / fUpper;
		worldFingers[i].y = (whiteLine[0].y - fingertips[i].y) * zUpper / fUpper;
	}
}

void frontToWorld(vector<Point> fingertips, Point vanishingPoint, vector<Point> whiteLine, vector<Point> blackLine, Mat image, Mat imageCut, vector<Point>& worldFingers) {

	vector<float> crossRatios = vector<float>(fingertips.size());
	Vec4i lineWhiteKey = Vec4i(whiteLine[0].x, whiteLine[0].y, whiteLine[1].x, whiteLine[1].y);
	Vec4i lineBlackKey = Vec4i(blackLine[0].x, blackLine[0].y, blackLine[1].x, blackLine[1].y);
	worldFingers = vector<Point>(fingertips.size());

	int counter = 0;
	for (int i = 0; i < fingertips.size(); i++) {
		Vec4i lineFingerTipVP = Vec4i(fingertips[i].x, fingertips[i].y, vanishingPoint.x, vanishingPoint.y);
		Point whiteIntersection = findIntersection(lineWhiteKey, lineFingerTipVP);
		Point blackIntersection = findIntersection(lineBlackKey, lineFingerTipVP);
		crossRatios[i] = (distance(blackIntersection, whiteIntersection) / distance(blackIntersection, vanishingPoint)) /
			(distance(fingertips[i], whiteIntersection) / distance(fingertips[i], vanishingPoint));
		counter++;
		worldFingers[i].y = 500 / crossRatios[i];
	}

}

Mat rotate(Mat src, double angle)
{
	Mat dst;
	Point2f pt(src.cols / 2., src.rows / 2.);
	Mat r = getRotationMatrix2D(pt, angle, 1.0);
	warpAffine(src, dst, r, Size(src.cols, src.rows));
	return dst;
}

void indexFingers(vector<Vec4i> topReference, vector<Vec4i> frontReference, vector<Point> fingertips, vector<Point> fingertips2, Mat imageCut, Mat imageCut2, vector<int>& topIndex, vector<int>& frontIndex, vector<float>& topKeyProportion, vector<float>& frontKeyProportion) {
	topIndex = vector<int>(fingertips.size());
	frontIndex = vector<int>(fingertips2.size());
	topKeyProportion = vector<float>(fingertips.size());
	frontKeyProportion = vector<float>(fingertips2.size());

	for (int i = 0; i < fingertips.size(); i++) {
		int closestX = camera_width;
		int index;
		for (int j = 0; j < topReference.size(); j++) {
			if (distancePointLine(fingertips[i], topReference[j]) < closestX) {
				topIndex[i] = topReference.size() - (1 + j);
				topKeyProportion[i] = 0.5;
				closestX = distancePointLine(fingertips[i], topReference[j]);
			}
		}
		if (topIndex[i] != 0 && topIndex[i] != topReference.size() - 1) {
			if (distancePointLine(fingertips[i], topReference[topReference.size() - topIndex[i]]) > distancePointLine(fingertips[i], topReference[topReference.size() - (2 + topIndex[i])])) topIndex[i]++;
			topKeyProportion[i] = distancePointLine(fingertips[i], topReference[topReference.size() - (1 + topIndex[i])]) / (distancePointLine(fingertips[i], topReference[topReference.size() - (1 + topIndex[i])]) + distancePointLine(fingertips[i], topReference[topReference.size() - topIndex[i]]));
		}
	}
	for (int i = 0; i < fingertips2.size(); i++) {
		int closestL = camera_width;
		int index;
		for (int j = 0; j < frontReference.size(); j++) {
			if (distancePointLine(fingertips2[i], frontReference[j]) < closestL) {
				frontIndex[i] = j;
				frontKeyProportion[i] = 0.5;
				closestL = distancePointLine(fingertips2[i], frontReference[j]);
			}
		}
		if (frontIndex[i] != 0 && frontIndex[i] != frontReference.size() - 1) {
			if (getIntersection(Point(frontReference[frontIndex[i]][0], frontReference[frontIndex[i]][1]), Point(frontReference[frontIndex[i]][2], frontReference[frontIndex[i]][3]), fingertips2[i], 0).x < fingertips2[i].x) frontIndex[i]++;
			frontKeyProportion[i] = distancePointLine(fingertips2[i], frontReference[frontIndex[i]]) / (distancePointLine(fingertips2[i], frontReference[frontIndex[i]]) + distancePointLine(fingertips2[i], frontReference[frontIndex[i] - 1]));
		}
		frontIndex[i] = frontIndex[i] + 2; //equalising factor
	}

}

void matchFingers(vector<Point> worldFingers, vector<Point> worldFingers2, vector<Point> fingertips, vector<Point> fingertips2, vector<int> topIndex, vector<int> frontIndex, Mat imageCut, Mat image, vector<Vec4i> topReference, vector<Point> whiteLine, int redY, vector<float> topKeyProportion, vector<float> frontKeyProportion, vector<int>& whiteKeyPresses, vector<int> lastWhiteKeyPresses) {

	int detectedKeyCount = 0;
	for (int i = 0; i < worldFingers.size(); i++) {
		int match = -1;
		for (int j = 0; j < worldFingers2.size(); j++) {
			if (topIndex[i] == frontIndex[j]) {
				if (match == -1) {
					if (abs(topKeyProportion[i] - frontKeyProportion[j]) < 0.25) match = j; //if new key press
					else if (abs(topKeyProportion[i] - frontKeyProportion[j]) < 0.35) { //if the key was being pressed in the last frame
						for (int j = 0; j < lastWhiteKeyPresses.size(); j++) {
							if (lastWhiteKeyPresses[j] == topIndex[i]) {
								match = j;
							}
						}
					}
				}
				else if (abs(worldFingers[i].y - worldFingers2[j].y) < abs(worldFingers[i].y - worldFingers2[match].y)) {
					if (abs(topKeyProportion[i] - frontKeyProportion[j]) < 0.25) match = j; //if new key press
					else if (abs(topKeyProportion[i] - frontKeyProportion[j]) < 0.35) { //if the key was being pressed in the last frame
						for (int j = 0; j < lastWhiteKeyPresses.size(); j++) {
							if (lastWhiteKeyPresses[j] == topIndex[i]) {
								match = j;
							}
						}
					}
				}
			}
		}

		//do match stuff
		if (match != -1) {
			int difference = worldFingers[i].y - worldFingers2[match].y;
			//putText(image, to_string(worldFingers[i].y), Point(fingertips[i].x, fingertips[i].y), FONT_HERSHEY_SIMPLEX, 1, Scalar(0, 0, 255), 2, LINE_AA);
			//putText(image, to_string(worldFingers2[match].y), Point(fingertips[i].x, fingertips[i].y - 40), FONT_HERSHEY_SIMPLEX, 1, Scalar(255, 0, 0), 2, LINE_AA);

			if (difference < -50) {
				rectangle(image, Point(topReference[topReference.size() - topIndex[i] - 1][0], redY),
					Point(topReference[topReference.size() - topIndex[i]][0], whiteLine[0].y),
					Scalar(168, 250, 50), 2);
				whiteKeyPresses[detectedKeyCount] = topIndex[i];
				detectedKeyCount++;
			}
			else if (difference < -25) { //if the key was being pressed in the last frame
				for (int j = 0; j < lastWhiteKeyPresses.size(); j++) {
					if (lastWhiteKeyPresses[j] == topIndex[i]) {
						rectangle(image, Point(topReference[topReference.size() - topIndex[i] - 1][0], redY),
							Point(topReference[topReference.size() - topIndex[i]][0], whiteLine[0].y),
							Scalar(168, 250, 50), 2);
						whiteKeyPresses[detectedKeyCount] = topIndex[i];
						detectedKeyCount++;
					}
				}
			}
		}
	}
	whiteKeyPresses.resize(detectedKeyCount);
}

void createBlackFilter(Mat& redLines, Mat init, vector<Point> blackLine, vector<Point> whiteLine, int& counter, int& redY, int& whiteKeyLength, vector<Point>& blackKeyLocations) {
	if (counter < 10) {
		Mat initYCrCb;
		Mat initHLS;
		Mat redHLS;
		Mat redYCrCb;
		Mat newLines;
		cvtColor(init, initYCrCb, COLOR_BGR2YCrCb);
		cvtColor(init, initHLS, COLOR_BGR2HLS);
		inRange(initYCrCb, Scalar(R_Y_MIN, R_Cr_MIN, R_Cb_MIN), Scalar(R_Y_MAX, R_Cr_MAX, R_Cb_MAX), redYCrCb);
		inRange(initHLS, Scalar(R_H_MIN, R_L_MIN, R_S_MIN), Scalar(R_H_MAX, R_L_MAX, R_S_MAX), redHLS);
		bitwise_and(redYCrCb, redHLS, newLines);
		erode(newLines, newLines, getStructuringElement(MORPH_RECT, Size(7, 7), Point(-1, -1)));
		dilate(newLines, newLines, getStructuringElement(MORPH_RECT, Size(4, 4), Point(-1, -1)));
		int redTop = blackLine[0].y - 2 * (whiteLine[0].y - blackLine[0].y);
		rectangle(newLines, Point(0, 0), Point(camera_width, redTop), 0, -1);
		rectangle(newLines, Point(0, blackLine[0].y), Point(camera_width, camera_height), 0, -1);
		if (counter == 0) newLines.copyTo(redLines);
		else bitwise_and(newLines, redLines, redLines);
		vector<vector<Point>> contours;
		findContours(redLines, contours, RETR_TREE, CHAIN_APPROX_SIMPLE);
		blackKeyLocations = vector<Point>(contours.size());
		int keyCounter = 0;
		for (int i = 0; i < contours.size(); i++) {
			if (contourArea(contours[i]) > 80) {
				Rect bound = boundingRect(contours[i]);
				blackKeyLocations[keyCounter] = (bound.br() + bound.tl()) / 2;
				if (keyCounter == 0) {
					redY = bound.y;
					whiteKeyLength = whiteLine[0].y - redY;
				}
				keyCounter++;
			}
		}
		blackKeyLocations.resize(keyCounter);
		sort(blackKeyLocations.begin(), blackKeyLocations.end(), sortHelper);
		counter++;
	}
}

void checkBlack(Mat image, Mat redLines, vector<Point> whiteLine, vector<Point> blackLine, vector<Point> fingertips, vector<Point> blackKeyLocations, vector<int>& blackKeyPresses) {
	Mat imageCopy;
	Mat YCrCb;
	Mat HLS;
	image.copyTo(imageCopy);
	Mat newHLS;
	Mat newYCrCb;
	Mat newRed;
	cvtColor(imageCopy, YCrCb, COLOR_BGR2YCrCb);
	cvtColor(imageCopy, HLS, COLOR_BGR2HLS);
	inRange(YCrCb, Scalar(R_Y_MIN, R_Cr_MIN, R_Cb_MIN), Scalar(R_Y_MAX, R_Cr_MAX, R_Cb_MAX), newYCrCb);
	inRange(HLS, Scalar(R_H_MIN, R_L_MIN, R_S_MIN), Scalar(R_H_MAX, R_L_MAX, R_S_MAX), newHLS);
	bitwise_and(newYCrCb, newHLS, newRed);
	erode(newRed, newRed, getStructuringElement(MORPH_RECT, Size(7, 7), Point(-1, -1)));
	dilate(newRed, newRed, getStructuringElement(MORPH_RECT, Size(4, 4), Point(-1, -1)));
	int redTop = blackLine[0].y - 2 * (whiteLine[0].y - blackLine[0].y);
	rectangle(newRed, Point(0, 0), Point(camera_width, redTop), 0, -1);
	rectangle(newRed, Point(0, blackLine[0].y), Point(camera_width, camera_height), 0, -1);
	Mat filtered;
	redLines.copyTo(filtered);
	filtered = filtered - newRed;
	vector<vector<Point>> contours;

	findContours(filtered, contours, RETR_TREE, CHAIN_APPROX_SIMPLE);
	blackKeyPresses = vector<int>(blackKeyLocations.size());
	int detectedKeyCount = 0;
	for (int i = 0; i < contours.size(); i++) {
		if (contourArea(contours[i]) > 80) {
			Rect bound = boundingRect(contours[i]);
			for (int z = 0; z < fingertips.size(); z++) { //changed i to z not tested
				if (fingertips[z].x > bound.x - 30 && fingertips[z].x < bound.x + bound.width + 30 &&
					fingertips[z].y > bound.y && fingertips[z].y < blackLine[0].y) {
					rectangle(image, Point(bound.x, bound.y), Point(bound.x + bound.width, blackLine[0].y), Scalar(0, 151, 252), 3);
					for (int j = 0; j < blackKeyLocations.size(); j++) {
						if (blackKeyLocations[j].x > bound.x && blackKeyLocations[j].x < bound.x + bound.width) {
							blackKeyPresses[detectedKeyCount] = j;
							detectedKeyCount++;
						}
					}
				}
			}
		}
	}
	blackKeyPresses.resize(detectedKeyCount);
}
/*
*    Title: Get pixel RGB value from webcam video in OpenCV
*    Author: Najam R. Syed
*    Date: 13.02.2018
*    Code version: 1.0
*    Availability: https://nrsyed.com/2018/02/12/get-pixel-rgb-value-from-webcam-video-in-opencv-c-and-python/
*/
void on_mouse_click(int event, int x, int y, int flags, void* ptr) {
	if (event == cv::EVENT_LBUTTONDOWN) {
		cv::Mat* snapshot = (cv::Mat*)ptr;
		cv::Vec3b pixel = snapshot->at<cv::Vec3b>(cv::Point(x, y));
		int b, g, r;
		b = pixel[0];
		g = pixel[1];
		r = pixel[2];
		std::string rgbText = "[" + std::to_string(r) + ", " + std::to_string(g)
			+ ", " + std::to_string(b) + "]";
		float luminance = 1 - (0.299 * r + 0.587 * g + 0.114 * b) / 255;
		cv::Scalar textColor;
		if (luminance < 0.5) {
			textColor = cv::Scalar(0, 0, 0);
		}
		else {
			textColor = cv::Scalar(255, 255, 255);
		}

		cv::Mat colorArray;
		colorArray = cv::Mat(80, 250, CV_8UC3, cv::Scalar(b, g, r));
		cv::putText(colorArray, rgbText, cv::Point2d(20, 80 - 20),
			cv::FONT_HERSHEY_SIMPLEX, 0.8, textColor);
		cv::imshow("Color", colorArray);
	}
}
/*
*    End of citation
*/
void createWhiteTest(Mat input, Mat& output) {
	input.copyTo(output);
	for (int i = 0; i < output.rows; i++) {
		for (int j = 0; j < output.cols; j++) {
			Vec3b& color = output.at<Vec3b>(i, j);
			if (abs(color[0] - color[1]) < 35 && abs(color[1] - color[2]) < 35 && abs(color[2] - color[0]) < 35) {
				color = Vec3b(0, 0, 0);
			}
		}
	}
}

void playNote(MidiTrack& trk, vector<int> newBlackKeyPresses, vector<int> removedBlackKeys, vector<int> newWhiteKeyPresses, vector<int> removedWhiteKeys, int counter, int& lastNote) {
	int elapse = (counter - lastNote) * 4;

	for (int i = 0; i < removedWhiteKeys.size(); i++) {
		switch (removedWhiteKeys[i]) {

		case 0: trk.noteOff(elapse, 0x53, 0x40);
			elapse = 0;
			break;
		case 1: trk.noteOff(elapse, 0x51, 0x40);
			elapse = 0;
			break;
		case 2: trk.noteOff(elapse, 0x4F, 0x40);
			elapse = 0;
			break;
		case 3: trk.noteOff(elapse, 0x4D, 0x40);
			elapse = 0;
			break;
		case 4: trk.noteOff(elapse, 0x4C, 0x40);
			elapse = 0;
			break;
		case 5: trk.noteOff(elapse, 0x4A, 0x40);
			elapse = 0;
			break;
		case 6: trk.noteOff(elapse, 0x48, 0x40);
			elapse = 0;
			break;
		case 7: trk.noteOff(elapse, 0x47, 0x40);
			elapse = 0;
			break;
		case 8: trk.noteOff(elapse, 0x45, 0x40);
			elapse = 0;
			break;
		case 9: trk.noteOff(elapse, 0x43, 0x40);
			elapse = 0;
			break;
		case 10:trk.noteOff(elapse, 0x41, 0x40);
			elapse = 0;
			break;
		case 11:trk.noteOff(elapse, 0x40, 0x40);
			elapse = 0;
			break;
		case 12:trk.noteOff(elapse, 0x3E, 0x40);
			elapse = 0;
			break;
		case 13:trk.noteOff(elapse, 0x3C, 0x40);
			elapse = 0;
			break;
		case 14:trk.noteOff(elapse, 0x3B, 0x40);
			elapse = 0;
			break;
		case 15:trk.noteOff(elapse, 0x39, 0x40);
			elapse = 0;
			break;
		case 16:trk.noteOff(elapse, 0x37, 0x40);
			elapse = 0;
			break;
		case 17:trk.noteOff(elapse, 0x35, 0x40);
			elapse = 0;
			break;
		case 18:trk.noteOff(elapse, 0x34, 0x40);
			elapse = 0;
			break;
		case 19:trk.noteOff(elapse, 0x32, 0x40);
			elapse = 0;
			break;
		case 20:trk.noteOff(elapse, 0x30, 0x40);
			elapse = 0;
			break;
		case 21:trk.noteOff(elapse, 0x2F, 0x40);
			elapse = 0;
			break;

		default: return;
		}
	}

	for (int i = 0; i < removedBlackKeys.size(); i++) {
		switch (removedBlackKeys[i]) {

		case 0: trk.noteOff(elapse, 0x2C, 0x40);
			elapse = 0;
			break;
		case 1: trk.noteOff(elapse, 0x2E, 0x40);
			elapse = 0;
			break;
		case 2: trk.noteOff(elapse, 0x31, 0x40);
			elapse = 0;
			break;
		case 3: trk.noteOff(elapse, 0x33, 0x40);
			elapse = 0;
			break;
		case 4: trk.noteOff(elapse, 0x36, 0x40);
			elapse = 0;
			break;
		case 5: trk.noteOff(elapse, 0x38, 0x40);
			elapse = 0;
			break;
		case 6: trk.noteOff(elapse, 0x3A, 0x40);
			elapse = 0;
			break;
		case 7: trk.noteOff(elapse, 0x3D, 0x40);
			elapse = 0;
			break;
		case 8: trk.noteOff(elapse, 0x3F, 0x40);
			elapse = 0;
			break;
		case 9: trk.noteOff(elapse, 0x42, 0x40);
			elapse = 0;
			break;
		case 10:trk.noteOff(elapse, 0x44, 0x40);
			elapse = 0;
			break;
		case 11:trk.noteOff(elapse, 0x46, 0x40);
			elapse = 0;
			break;
		case 12:trk.noteOff(elapse, 0x49, 0x40);
			elapse = 0;
			break;
		case 13:trk.noteOff(elapse, 0x4B, 0x40);
			elapse = 0;
			break;
		case 14:trk.noteOff(elapse, 0x4E, 0x40);
			elapse = 0;
			break;
		case 15:trk.noteOff(elapse, 0x50, 0x40);
			elapse = 0;
			break;
		case 16:trk.noteOff(elapse, 0x52, 0x40);
			elapse = 0;
			break;

		default: return;
		}
	}

	for (int i = 0; i < newWhiteKeyPresses.size(); i++) {
		switch (newWhiteKeyPresses[i]) {

		case 0: trk.noteOn(elapse, 0x53, 0x40);
			elapse = 0;
			break;
		case 1: trk.noteOn(elapse, 0x51, 0x40);
			elapse = 0;
			break;
		case 2: trk.noteOn(elapse, 0x4F, 0x40);
			elapse = 0;
			break;
		case 3: trk.noteOn(elapse, 0x4D, 0x40);
			elapse = 0;
			break;
		case 4: trk.noteOn(elapse, 0x4C, 0x40);
			elapse = 0;
			break;
		case 5: trk.noteOn(elapse, 0x4A, 0x40);
			elapse = 0;
			break;
		case 6: trk.noteOn(elapse, 0x48, 0x40);
			elapse = 0;
			break;
		case 7: trk.noteOn(elapse, 0x47, 0x40);
			elapse = 0;
			break;
		case 8: trk.noteOn(elapse, 0x45, 0x40);
			elapse = 0;
			break;
		case 9: trk.noteOn(elapse, 0x43, 0x40);
			elapse = 0;
			break;
		case 10:trk.noteOn(elapse, 0x41, 0x40);
			elapse = 0;
			break;
		case 11:trk.noteOn(elapse, 0x40, 0x40);
			elapse = 0;
			break;
		case 12:trk.noteOn(elapse, 0x3E, 0x40);
			elapse = 0;
			break;
		case 13:trk.noteOn(elapse, 0x3C, 0x40);
			elapse = 0;
			break;
		case 14:trk.noteOn(elapse, 0x3B, 0x40);
			elapse = 0;
			break;
		case 15:trk.noteOn(elapse, 0x39, 0x40);
			elapse = 0;
			break;
		case 16:trk.noteOn(elapse, 0x37, 0x40);
			elapse = 0;
			break;
		case 17:trk.noteOn(elapse, 0x35, 0x40);
			elapse = 0;
			break;
		case 18:trk.noteOn(elapse, 0x34, 0x40);
			elapse = 0;
			break;
		case 19:trk.noteOn(elapse, 0x32, 0x40);
			elapse = 0;
			break;
		case 20:trk.noteOn(elapse, 0x30, 0x40);
			elapse = 0;
			break;
		case 21:trk.noteOn(elapse, 0x2F, 0x40);
			elapse = 0;
			break;

		default: return;
		}
	}

	for (int i = 0; i < newBlackKeyPresses.size(); i++) {
		switch (newBlackKeyPresses[i]) {

		case 0: trk.noteOn(elapse, 0x2C, 0x40);
			elapse = 0;
			break;
		case 1: trk.noteOn(elapse, 0x2E, 0x40);
			elapse = 0;
			break;
		case 2: trk.noteOn(elapse, 0x31, 0x40);
			elapse = 0;
			break;
		case 3: trk.noteOn(elapse, 0x33, 0x40);
			elapse = 0;
			break;
		case 4: trk.noteOn(elapse, 0x36, 0x40);
			elapse = 0;
			break;
		case 5: trk.noteOn(elapse, 0x38, 0x40);
			elapse = 0;
			break;
		case 6: trk.noteOn(elapse, 0x3A, 0x40);
			elapse = 0;
			break;
		case 7: trk.noteOn(elapse, 0x3D, 0x40);
			elapse = 0;
			break;
		case 8: trk.noteOn(elapse, 0x3F, 0x40);
			elapse = 0;
			break;
		case 9: trk.noteOn(elapse, 0x42, 0x40);
			elapse = 0;
			break;
		case 10:trk.noteOn(elapse, 0x44, 0x40);
			elapse = 0;
			break;
		case 11:trk.noteOn(elapse, 0x46, 0x40);
			elapse = 0;
			break;
		case 12:trk.noteOn(elapse, 0x49, 0x40);
			elapse = 0;
			break;
		case 13:trk.noteOn(elapse, 0x4B, 0x40);
			elapse = 0;
			break;
		case 14:trk.noteOn(elapse, 0x4E, 0x40);
			elapse = 0;
			break;
		case 15:trk.noteOn(elapse, 0x50, 0x40);
			elapse = 0;
			break;
		case 16:trk.noteOn(elapse, 0x52, 0x40);
			elapse = 0;
			break;

		default: return;
		}
	}
	if (elapse == 0) lastNote = counter;
}

void checkNoteChanges(vector<int>& blackKeyPresses, vector<int>& whiteKeyPresses, vector<int>& lastBlackKeyPresses, vector<int>& lastWhiteKeyPresses, int counter, int& lastNote, MidiTrack& trk) {
	vector<int> newBlackKeyPresses = vector<int>(blackKeyPresses);
	int ct = 0;
	for (int i = 0; i < blackKeyPresses.size(); i++) {
		bool noMatch = 1;
		for (int j = 0; j < lastBlackKeyPresses.size(); j++) {
			if (blackKeyPresses[i] == lastBlackKeyPresses[j]) noMatch = 0;
		}
		if (noMatch) {
			newBlackKeyPresses[ct] = blackKeyPresses[i];
			ct++;
		}
	}
	newBlackKeyPresses.resize(ct);

	vector<int> removedBlackKeys = vector<int>(lastBlackKeyPresses);
	ct = 0;
	for (int i = 0; i < lastBlackKeyPresses.size(); i++) {
		bool noMatch = 1;
		for (int j = 0; j < blackKeyPresses.size(); j++) {
			if (lastBlackKeyPresses[i] == blackKeyPresses[j]) noMatch = 0;
		}
		if (noMatch) {
			removedBlackKeys[ct] = lastBlackKeyPresses[i];
			ct++;
		}
	}
	removedBlackKeys.resize(ct);

	vector<int> newWhiteKeyPresses = vector<int>(whiteKeyPresses);
	ct = 0;
	for (int i = 0; i < whiteKeyPresses.size(); i++) {
		bool noMatch = 1;
		for (int j = 0; j < lastWhiteKeyPresses.size(); j++) {
			if (whiteKeyPresses[i] == lastWhiteKeyPresses[j]) noMatch = 0;
		}
		if (noMatch) {
			newWhiteKeyPresses[ct] = whiteKeyPresses[i];
			ct++;
		}
	}
	newWhiteKeyPresses.resize(ct);

	vector<int> removedWhiteKeys = vector<int>(lastWhiteKeyPresses);
	ct = 0;
	for (int i = 0; i < lastWhiteKeyPresses.size(); i++) {
		bool noMatch = 1;
		for (int j = 0; j < whiteKeyPresses.size(); j++) {
			if (lastWhiteKeyPresses[i] == whiteKeyPresses[j]) noMatch = 0;
		}
		if (noMatch) {
			removedWhiteKeys[ct] = lastWhiteKeyPresses[i];
			ct++;
		}
	}
	removedWhiteKeys.resize(ct);

	playNote(trk, newBlackKeyPresses, removedBlackKeys, newWhiteKeyPresses, removedWhiteKeys, counter, lastNote);

	lastBlackKeyPresses = blackKeyPresses;
	lastWhiteKeyPresses = whiteKeyPresses;
}

void closeMidi(MidiHeader hd, MidiTrack trk) {
	trk.endOfTrack(0);
	MidiFile output("HAPPY_FARMER2.mid", hd, trk);
	output.closeFile();
	Sleep(5000);
}



int main() {
	VideoCapture cap;
	VideoCapture cap2;

	Mat init, initYCrCb, keyYCrCb, keyBW, cannyKey, keyCut, redLines;
	Mat init2, initYCrCb2, initHLS2, keyYCrCb2, keyHLS2, keyFiltered2, keyBW2, cannyKey2, keyCut2;

	Point vanishingPoint;
	vector<Point> fingertips;
	vector<Point> worldFingers;
	vector<Point> fingertips2;
	vector<Point> worldFingers2;

	vector<Point> mFingertips;
	vector<Point> mWorldFingers;
	vector<Point> mFingertips2;
	vector<Point> mWorldFingers2;

	vector<Point> whiteLine = vector<Point>(2);
	vector<Point> blackLine = vector<Point>(2);
	vector<Point> whiteLine2 = vector<Point>(2);
	vector<Point> blackLine2 = vector<Point>(2);

	vector<Point> blackKeyLocations;
	int redY;

	//meta events

	uint16_t division = 4;
	unsigned int quarterNote = 20;
	unsigned int halfNote = 40;
	MidiHeader hd(0, 1, division);
	MidiTrack trk;

	trk.setTempo(0, 50000); //120bpm
	trk.programChange(0, 1); //sets instrument as piano
	trk.keySignature(0, 0, 0); //C Major

	cap.open("trial 2 upper.mp4");
	cap2.open("trial 2 front.mp4");
	namedWindow("final", WINDOW_NORMAL);
	namedWindow("front", WINDOW_NORMAL);
	resizeWindow("final", camera_width * 2 / 3, camera_height * 2 / 3);
	resizeWindow("front", camera_width * 2 / 3, camera_height * 2 / 3);

	createTrackbars();

	cap >> init;
	cvtColor(init, initYCrCb, COLOR_BGR2YCrCb);
	inRange(initYCrCb, Scalar(K_Y_MIN, K_Cr_MIN, K_Cb_MIN), Scalar(K_Y_MAX, K_Cr_MAX, K_Cb_MAX), keyYCrCb);
	bitwise_and(init, init, keyCut, keyYCrCb);
	cvtColor(keyCut, keyBW, COLOR_BGR2GRAY);
	Canny(keyBW, cannyKey, 50, 200, 3);
	Hough(cannyKey, vanishingPoint, whiteLine, blackLine, 0);

	rectangle(init, Point(0, 0), Point(camera_width, blackLine[0].y), 0, -1);
	rectangle(init, Point(0, whiteLine[0].y), Point(camera_width, camera_height), 0, -1);
	cvtColor(init, initYCrCb, COLOR_BGR2YCrCb);
	inRange(initYCrCb, Scalar(K_Y_MIN, K_Cr_MIN, K_Cb_MIN), Scalar(K_Y_MAX, K_Cr_MAX, K_Cb_MAX), keyYCrCb);
	Mat newCut;
	bitwise_and(init, init, newCut, keyYCrCb);
	cvtColor(newCut, keyBW, COLOR_BGR2GRAY);
	Canny(keyBW, cannyKey, 50, 200, 3);


	Mat keyLine(init.rows, init.cols, 1, Scalar(0));
	line(keyLine, Point(0, (blackLine[0].y + whiteLine[0].y) / 2), Point(camera_width, (blackLine[0].y + whiteLine[0].y) / 2), Scalar(256, 256, 256), 1, LINE_AA);

	Mat topKeyLocations;
	bitwise_and(keyLine, keyLine, topKeyLocations, cannyKey);
	cvtColor(cannyKey, cannyKey, COLOR_GRAY2BGR);
	line(cannyKey, Point(0, (blackLine[0].y + whiteLine[0].y) / 2), Point(camera_width, (blackLine[0].y + whiteLine[0].y) / 2), Scalar(0, 0, 256), 1, LINE_AA);

	vector<Point> topKeyPoints;
	findNonZero(topKeyLocations, topKeyPoints);

	sort(topKeyPoints.begin(), topKeyPoints.end(), sortHelper);

	vector<Vec4i> topReference(topKeyPoints.size());
	int ct = 0;
	int lastX = 0;
	for (int i = 0; i < topKeyPoints.size(); i++) {
		if (topKeyPoints[i].x - lastX > camera_width / 30) {
			line(cannyKey, Point(topKeyPoints[i].x, 0), Point(topKeyPoints[i].x, camera_height), Scalar(256, 0, 0), 1, LINE_AA);
			topReference[ct] = Vec4i(topKeyPoints[i].x, 0, topKeyPoints[i].x, camera_height);
			lastX = topKeyPoints[i].x;
			ct++;
		}
	}
	if (ct > 0) topReference.resize(ct - 1);

	cap2 >> init2;
	cvtColor(init2, initYCrCb2, COLOR_BGR2YCrCb);
	cvtColor(init2, initHLS2, COLOR_BGR2HLS);
	inRange(initYCrCb2, Scalar(K_Y2_MIN, K_Cr2_MIN, K_Cb2_MIN), Scalar(K_Y2_MAX, K_Cr2_MAX, K_Cb2_MAX), keyYCrCb2);
	inRange(initHLS2, Scalar(K_H_MIN, K_L_MIN, K_S_MIN), Scalar(K_H_MAX, K_L_MAX, K_S_MAX), keyHLS2);
	bitwise_and(keyHLS2, keyYCrCb2, keyFiltered2);
	erode(keyFiltered2, keyFiltered2, getStructuringElement(MORPH_RECT, Size(7, 7), Point(-1, -1)));
	bitwise_and(init2, init2, keyCut2, keyFiltered2);
	cvtColor(keyCut2, keyBW2, COLOR_BGR2GRAY);
	Canny(keyBW2, cannyKey2, 50, 200, 3);
	Hough(cannyKey2, vanishingPoint, whiteLine2, blackLine2, 1);

	rectangle(init2, Point(0, 0), Point(camera_width, whiteLine2[0].y), 0, -1);
	rectangle(init2, Point(0, blackLine2[0].y), Point(camera_width, camera_height), 0, -1);
	cvtColor(init2, initYCrCb2, COLOR_BGR2YCrCb);
	cvtColor(init2, initHLS2, COLOR_BGR2HLS);
	inRange(initYCrCb2, Scalar(K_Y2_MIN, K_Cr2_MIN, K_Cb2_MIN), Scalar(K_Y2_MAX, K_Cr2_MAX, K_Cb2_MAX), keyYCrCb2);
	inRange(initHLS2, Scalar(K_H_MIN, K_L_MIN, K_S_MIN), Scalar(K_H_MAX, K_L_MAX, K_S_MAX), keyHLS2);
	bitwise_and(keyHLS2, keyYCrCb2, keyFiltered2);
	erode(keyFiltered2, keyFiltered2, getStructuringElement(MORPH_RECT, Size(7, 7), Point(-1, -1)));
	Mat newCut2;
	bitwise_and(init2, init2, newCut2, keyFiltered2);
	cvtColor(newCut2, keyBW2, COLOR_BGR2GRAY);
	Canny(keyBW2, cannyKey2, 50, 200, 3);
	vector<Vec4i> linesP2;
	HoughLinesP(cannyKey2, linesP2, 1, CV_PI / 180, 40, 40, 10);
	sort(linesP2.begin(), linesP2.end(), sortHelper2);

	vector<Vec4i> frontReference(linesP2.size());
	ct = 0;
	lastX = 0;
	cvtColor(cannyKey2, cannyKey2, COLOR_GRAY2BGR);
	for (size_t i = 0; i < linesP2.size(); i++) {
		Vec4i l = linesP2[i];
		Point first = Point(l[0], l[1]);
		Point second = Point(l[2], l[3]);
		float sl = slope(first, second);
		Point start, finish;
		if (l[0] - lastX > camera_width / 50 && abs(sl) > 1) {
			if (sl == 100000000000) {
				start.x = l[0];
				finish.x = l[0];
				start.y = 0;
				finish.y = camera_height;
			}
			else {
				start.x = first.x - first.y / sl;
				finish.x = -(second.y - camera_height) / sl + second.x;
				start.y = 0;
				finish.y = camera_height;
			}
			line(cannyKey2, start, finish, Scalar(256, 0, 0), 1, LINE_AA);
			lastX = l[0];
			frontReference[ct] = Vec4i(start.x, start.y, finish.x, finish.y);
			ct++;
		}
	}
	frontReference.resize(ct);

	int filtercount = 0;
	int whiteKeyLength;
	int framecounter = 0;
	int lastNote = 0;
	vector<int> lastWhiteKeyPresses;
	vector<int> lastBlackKeyPresses;
	while (1) {

		framecounter++;
		Mat image;
		Mat imageYCrCb;
		Mat imageFiltered;

		Mat image2;
		Mat imageYCrCb2;
		Mat imageHLS;
		Mat imageFiltered2;
		Mat YCrCbFilter;
		Mat HLSFilter;
		Mat BGRFilter;
		Mat imageCut2, imageCut;
		Mat whiteDiff;

		cap >> image;

		cap2 >> image2;
		if (image.empty() || image2.empty()) break;
		createBlackFilter(redLines, image, blackLine, whiteLine, filtercount, redY, whiteKeyLength, blackKeyLocations);
		//color filters
		rectangle(image2, Point(0, 0), Point(camera_width, camera_height / 2 + 20), 0, -1);
		cvtColor(image, imageYCrCb, COLOR_BGR2YCrCb);
		inRange(imageYCrCb, Scalar(Y_MIN, Cr_MIN, Cb_MIN), Scalar(Y_MAX, Cr_MAX, Cb_MAX), imageFiltered);
		erode(imageFiltered, imageFiltered, getStructuringElement(MORPH_RECT, Size(7, 7), Point(-1, -1)));

		cvtColor(image2, imageYCrCb2, COLOR_BGR2YCrCb);
		cvtColor(image2, imageHLS, COLOR_BGR2HLS);
		inRange(imageYCrCb2, Scalar(Y2_MIN, Cr2_MIN, Cb2_MIN), Scalar(Y2_MAX, Cr2_MAX, Cb2_MAX), YCrCbFilter);
		inRange(imageHLS, Scalar(H_MIN, L_MIN, S_MIN), Scalar(H_MAX, L_MAX, S_MAX), HLSFilter);
		createWhiteTest(image2, whiteDiff);
		inRange(whiteDiff, Scalar(B_MIN, G_MIN, R_MIN), Scalar(B_MAX, G_MAX, R_MAX), BGRFilter);
		bitwise_and(HLSFilter, YCrCbFilter, imageFiltered2);
		bitwise_and(BGRFilter, imageFiltered2, imageFiltered2);
		erode(imageFiltered2, imageFiltered2, getStructuringElement(MORPH_RECT, Size(7, 7), Point(-1, -1)));

		//result of color filters, application of canny edges on said result
		bitwise_and(image, image, imageCut, imageFiltered);
		bitwise_and(image2, image2, imageCut2, imageFiltered2);

		//calculation of fingertips and display on the result of color filter result
		calculateFingers(imageFiltered, imageCut, fingertips, 1, whiteLine); //convexity defects method
		calculateFingers(imageFiltered2, imageCut2, fingertips2, -1, whiteLine2);

		calculateFingers2(imageFiltered, imageCut, mFingertips, 1, whiteLine); //local extrema method
		calculateFingers2(imageFiltered2, imageCut2, mFingertips2, -1, whiteLine2);

		vector<int> topIndex;
		vector<int> frontIndex;
		vector<int> topIndex2;
		vector<int> frontIndex2;
		vector<float> topKeyProportion;
		vector<float> frontKeyProportion;
		vector<float> topKeyProportion2;
		vector<float> frontKeyProportion2;
		//index fingertips according to piano keys
		indexFingers(topReference, frontReference, fingertips, fingertips2, imageCut, imageCut2, topIndex, frontIndex, topKeyProportion, frontKeyProportion);
		indexFingers(topReference, frontReference, mFingertips, mFingertips2, imageCut, imageCut2, topIndex2, frontIndex2, topKeyProportion2, frontKeyProportion2);
		//combine fingertips from convexity defects and local extrema methods
		fillMissingFingers(image, imageCut, imageCut2, fingertips, fingertips2, mFingertips, mFingertips2, topIndex, frontIndex, topIndex2, frontIndex2, topKeyProportion, frontKeyProportion, topKeyProportion2, frontKeyProportion2, image2);
		//calculate world coordinates
		upperToWorld(fingertips, imageCut, whiteLine, worldFingers, whiteKeyLength);
		frontToWorld(fingertips2, vanishingPoint, whiteLine2, blackLine2, image2, imageCut2, worldFingers2);
		vector<int> blackKeyPresses = vector<int>(10);
		//calculate pressed black keys
		checkBlack(image, redLines, whiteLine, blackLine, fingertips, blackKeyLocations, blackKeyPresses);
		vector<int> whiteKeyPresses = vector<int>(10);
		//match fingertips based on indexes and calculate pressed white keys
		matchFingers(worldFingers, worldFingers2, fingertips, fingertips2, topIndex, frontIndex, imageCut, image, topReference, whiteLine, redY, topKeyProportion, frontKeyProportion, whiteKeyPresses, lastWhiteKeyPresses);
		//create midi file based on the pressed piano keys
		checkNoteChanges(blackKeyPresses, whiteKeyPresses, lastBlackKeyPresses, lastWhiteKeyPresses, framecounter, lastNote, trk);
		//display fingertips and pressed keys for upper camera
		imshow("final", image);
		//display fingertips for front camera
		imshow("front", image2);

		waitKey(1);
	}

	closeMidi(hd, trk);
}