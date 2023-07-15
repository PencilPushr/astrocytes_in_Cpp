#ifndef REACTIONS_H
#define REACTIONS_H

class Reactions
{
public:

    Reaction(const std::vector<float>& ini, int number_of_rec, int number_of_var) :
        ini{ ini }, number_of_rec{ number_of_rec }, number_of_var{ number_of_var }
    {
    }

public:
    float reaction1(int i)
    {
        return rec(i, 0) * rec(i, 1);
    }

    float reaction2(int i)
    {
        return rec(i);
    }

    float reaction3(int i) {
        float rec7 = rec(i, 7);
        float rec14 = rec(i, 14);
        float rec5 = rec(i, 5);
        float rec4 = rec(i, 4);
        float rec15 = rec(i, 15);
        float rec17 = rec(i, 17);
        float rec26 = rec(i, 26);
        float rec9 = rec(i, 9);

        float numerator = std::pow(rec14, rec5) * std::pow(rec4, rec5);
        float denominator = (std::pow(rec4, rec5) + std::pow(rec14, rec5)) * (std::pow(rec4, rec5) + std::pow(rec15, rec5));
        float fraction1 = numerator / denominator;
        float fraction2 = rec17 / (rec17 + std::pow(rec16, rec26));
        float difference = rec9 - rec4;

        return 4 * rec7 * fraction1 * fraction2 * difference;
    }

private:
    // struct {} sigma1;
    // struct {} sigma2;
    std::vector<float>& ini;
    int number_of_rec;
    int number_of_var;
    float rec(int i, int x)
    {
        return ini[(i / numberofrec) * numberofvar + x];
    }
}

#endif // REACTIONS_H