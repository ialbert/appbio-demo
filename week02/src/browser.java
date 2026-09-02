import java.awt.*;
import java.awt.event.*;
import java.io.*;
import java.util.*;
import java.util.List;
import java.util.zip.GZIPInputStream;
import javax.swing.*;

/**
 * Minimal genome browser: sequence track + GFF track, drag to pan, wheel to zoom.
 */
public class browser extends JFrame {

    static final Color A = new Color(50, 180, 70);
    static final Color T = new Color(220, 55, 55);
    static final Color G = new Color(230, 175, 20);
    static final Color C = new Color(50, 90, 210);
    static final Color N = new Color(160, 160, 160);
    static final Color GENE = new Color(70, 130, 190);
    static final Color CDS = new Color(40, 75, 150);
    static final Color UTR = new Color(130, 175, 130);
    static final Color OTHER = new Color(150, 150, 160);
    static final Color BG = new Color(248, 248, 250);
    static final Color AXIS = new Color(50, 50, 55);

    public static void main(String[] args) {
        String fa = args.length > 0 ? args[0] : "fasta/Bundibugyo.fa";
        String gff = args.length > 1 ? args[1] : "gff/Bundibugyo.gff.gz";
        Fasta seq = Fasta.load(fa);
        List<Feat> feats = Gff.load(gff);
        SwingUtilities.invokeLater(new Runnable() {
            public void run() {
                new browser(seq, feats).setVisible(true);
            }
        });
    }

    static void die(String msg, Exception e) {
        System.err.println("error: " + msg);
        if (e != null) e.printStackTrace();
        System.exit(1);
    }

    browser(Fasta seq, List<Feat> feats) {
        super(seq.id + "  (" + seq.len() + " bp)");
        setDefaultCloseOperation(EXIT_ON_CLOSE);
        View v = new View(seq.len());
        Panel p = new Panel(seq, feats, v);
        JLabel status = new JLabel(" drag to pan   wheel to zoom   Fit resets");
        status.setBorder(BorderFactory.createEmptyBorder(4, 8, 4, 8));
        p.status = status;

        JButton fit = new JButton("Fit");
        JButton zin = new JButton("+");
        JButton zout = new JButton("-");
        fit.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                v.fit(p.trackW());
                p.repaint();
            }
        });
        zin.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                v.zoom(p.trackW() / 2, 1.25);
                v.clamp(p.trackW());
                p.repaint();
            }
        });
        zout.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                v.zoom(p.trackW() / 2, 0.8);
                v.clamp(p.trackW());
                p.repaint();
            }
        });

        JPanel bar = new JPanel(new FlowLayout(FlowLayout.LEFT, 6, 4));
        bar.add(fit);
        bar.add(zout);
        bar.add(zin);
        bar.add(new JLabel(seq.id + "   " + seq.len() + " bp"));

        add(bar, BorderLayout.NORTH);
        add(p, BorderLayout.CENTER);
        add(status, BorderLayout.SOUTH);
        setSize(1100, 400);
        setLocationRelativeTo(null);
    }
}

class Fasta {
    String id = "";
    String seq = "";

    int len() { return seq.length(); }

    char at(int i) {
        return (i < 0 || i >= seq.length()) ? 'N' : seq.charAt(i);
    }

    static Fasta load(String path) {
        Fasta f = new Fasta();
        StringBuilder b = new StringBuilder();
        BufferedReader r = Io.open(path);
        try {
            String line;
            while ((line = r.readLine()) != null) {
                if (line.startsWith(">")) {
                    if (f.id.isEmpty()) {
                        f.id = line.substring(1).trim().split("\\s+")[0];
                    }
                } else {
                    b.append(line.trim().toUpperCase());
                }
            }
            r.close();
        } catch (IOException e) {
            browser.die("reading FASTA " + path, e);
        }
        f.seq = b.toString();
        if (f.seq.isEmpty()) browser.die("empty FASTA: " + path, null);
        return f;
    }
}

class Feat {
    String type, name;
    int start, end; // 1-based inclusive
    char strand;

    Feat(String type, int start, int end, char strand, String name) {
        this.type = type;
        this.start = start;
        this.end = end;
        this.strand = strand;
        this.name = name;
    }

    boolean overlaps(int a, int b) {
        return start <= b && end >= a;
    }
}

class Gff {
    static final Set<String> SKIP = new HashSet<String>(Arrays.asList(
            "region", "exon", "mRNA", "regulatory_region",
            "polyA_signal_sequence", "sequence_feature"));

    static List<Feat> load(String path) {
        List<Feat> out = new ArrayList<Feat>();
        if (!new File(path).exists()) return out;
        BufferedReader r = Io.open(path);
        try {
            String line;
            while ((line = r.readLine()) != null) {
                if (line.isEmpty() || line.startsWith("#")) continue;
                String[] c = line.split("\t", -1);
                if (c.length < 8) continue;
                if (SKIP.contains(c[2])) continue;
                int start = Integer.parseInt(c[3]);
                int end = Integer.parseInt(c[4]);
                char strand = c[6].isEmpty() ? '.' : c[6].charAt(0);
                out.add(new Feat(c[2], start, end, strand, label(c)));
            }
            r.close();
        } catch (IOException e) {
            browser.die("reading GFF " + path, e);
        }
        Collections.sort(out, new Comparator<Feat>() {
            public int compare(Feat a, Feat b) {
                int d = a.start - b.start;
                return d != 0 ? d : a.end - b.end;
            }
        });
        return out;
    }

    static String label(String[] c) {
        String attrs = c.length > 8 ? c[8] : "";
        String name = attr(attrs, "gene");
        if (name == null) name = attr(attrs, "Name");
        if (name == null) name = attr(attrs, "product");
        if (name == null) name = c[2];
        return name;
    }

    static String attr(String attrs, String key) {
        for (String p : attrs.split(";")) {
            int eq = p.indexOf('=');
            if (eq < 0) continue;
            if (p.substring(0, eq).trim().equals(key)) {
                try {
                    return java.net.URLDecoder.decode(p.substring(eq + 1), "UTF-8");
                } catch (Exception e) {
                    return p.substring(eq + 1);
                }
            }
        }
        return null;
    }
}

class View {
    final int len;
    double start; // 0-based genomic coord at left edge of track
    double px = 1; // pixels per base
    static final double MAX_PX = 24;

    View(int len) {
        this.len = len;
        this.start = 0;
    }

    void fit(int w) {
        if (w <= 0) return;
        px = Math.max(0.01, (double) w / len);
        start = 0;
        clamp(w);
    }

    void pan(double dx, int w) {
        start -= dx / px;
        clamp(w);
    }

    void zoom(int mx, double factor) {
        double g = start + mx / px;
        px = Math.min(MAX_PX, Math.max(0.01, px * factor));
        start = g - mx / px;
    }

    void clamp(int w) {
        if (w <= 0) return;
        double vis = w / px;
        if (vis >= len) {
            start = 0;
            px = (double) w / len;
        } else {
            if (start < 0) start = 0;
            if (start + vis > len) start = len - vis;
        }
    }

    int xOf(double pos, int left) {
        return left + (int) Math.round((pos - start) * px);
    }

    double posOf(int x, int left) {
        return start + (x - left) / px;
    }

    int visStart() {
        return Math.max(1, (int) Math.floor(start) + 1);
    }

    int visEnd(int w) {
        return Math.min(len, (int) Math.ceil(start + w / px));
    }
}

class Panel extends JPanel implements MouseListener, MouseMotionListener, MouseWheelListener {
    final Fasta fa;
    final List<Feat> feats;
    final View v;
    JLabel status;
    int lastX;
    boolean fitted;
    Feat hover;

    static final int LEFT = 44;
    static final int AXIS = 28;
    static final int SEQ = 40;
    static final int GAP = 14;
    static final int LANE = 26;

    Panel(Fasta fa, List<Feat> feats, View v) {
        this.fa = fa;
        this.feats = feats;
        this.v = v;
        setBackground(browser.BG);
        setFocusable(true);
        addMouseListener(this);
        addMouseMotionListener(this);
        addMouseWheelListener(this);
        addComponentListener(new ComponentAdapter() {
            public void componentResized(ComponentEvent e) {
                if (!fitted && trackW() > 0) {
                    v.fit(trackW());
                    fitted = true;
                } else {
                    v.clamp(trackW());
                }
                repaint();
            }
        });
    }

    int trackW() {
        return Math.max(1, getWidth() - LEFT);
    }

    protected void paintComponent(Graphics g0) {
        super.paintComponent(g0);
        Graphics2D g = (Graphics2D) g0;
        g.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
        g.setRenderingHint(RenderingHints.KEY_TEXT_ANTIALIASING, RenderingHints.VALUE_TEXT_ANTIALIAS_ON);
        int w = trackW();
        v.clamp(w);
        drawAxis(g, w);
        drawSeq(g, w);
        drawGff(g, w);
    }

    void drawAxis(Graphics2D g, int w) {
        g.setColor(browser.AXIS);
        g.setFont(new Font("SansSerif", Font.PLAIN, 11));
        int y = 22;
        g.drawLine(LEFT, y, LEFT + w, y);
        int a = v.visStart();
        int b = v.visEnd(w);
        int step = tick(b - a);
        int first = ((a + step - 1) / step) * step;
        for (int p = first; p <= b; p += step) {
            int x = v.xOf(p - 1, LEFT);
            g.drawLine(x, y - 4, x, y + 4);
            g.drawString(fmt(p), x + 3, y - 6);
        }
    }

    void drawSeq(Graphics2D g, int w) {
        int y = AXIS;
        g.setColor(Color.WHITE);
        g.fillRect(LEFT, y, w, SEQ);
        g.setColor(new Color(220, 220, 224));
        g.drawLine(LEFT, y, LEFT + w, y);
        g.drawLine(LEFT, y + SEQ, LEFT + w, y + SEQ);
        label(g, "seq", y + SEQ / 2 + 4);

        boolean letters = v.px >= 9;
        int gap = v.px >= 10 ? 1 : 0;
        Shape old = g.getClip();
        g.clipRect(LEFT, y, w, SEQ);
        if (v.px >= 1) {
            int i0 = Math.max(0, (int) Math.floor(v.start));
            int i1 = Math.min(fa.len(), (int) Math.ceil(v.start + w / v.px));
            int bw = Math.max(1, (int) Math.round(v.px) - gap);
            Font font = new Font(Font.MONOSPACED, Font.BOLD, Math.min(16, Math.max(8, (int) v.px - 2)));
            g.setFont(font);
            FontMetrics fm = g.getFontMetrics();
            for (int i = i0; i < i1; i++) {
                char b = fa.at(i);
                int x = v.xOf(i, LEFT);
                g.setColor(baseColor(b));
                g.fillRect(x, y + 4, bw, SEQ - 8);
                if (letters) {
                    g.setColor(b == 'G' ? Color.DARK_GRAY : Color.WHITE);
                    String s = String.valueOf(b);
                    g.drawString(s, x + (bw - fm.stringWidth(s)) / 2,
                            y + SEQ / 2 + fm.getAscent() / 2 - 2);
                }
            }
        } else {
            for (int x = LEFT; x < LEFT + w; x++) {
                int i = (int) v.posOf(x, LEFT);
                g.setColor(baseColor(fa.at(i)));
                g.drawLine(x, y + 4, x, y + SEQ - 5);
            }
        }
        g.setClip(old);
    }

    void drawGff(Graphics2D g, int w) {
        int y0 = AXIS + SEQ + GAP;
        label(g, "gff", y0 + 14);

        int a = v.visStart();
        int b = v.visEnd(w);
        List<Feat> genes = new ArrayList<Feat>();
        List<Feat> cds = new ArrayList<Feat>();
        List<Feat> rest = new ArrayList<Feat>();
        for (Feat f : feats) {
            if (!f.overlaps(a, b)) continue;
            if (f.type.equals("gene")) genes.add(f);
            else if (f.type.equals("CDS")) cds.add(f);
            else rest.add(f);
        }

        Shape old = g.getClip();
        g.clipRect(LEFT, y0 - 2, w, 3 * LANE);
        drawLane(g, genes, y0, 22, true);
        drawLane(g, cds, y0 + LANE, 16, false);
        if (!rest.isEmpty()) drawLane(g, rest, y0 + 2 * LANE, 12, true);
        g.setClip(old);
    }

    void drawLane(Graphics2D g, List<Feat> row, int y, int h, boolean label) {
        Font font = new Font("SansSerif", Font.BOLD, 11);
        g.setFont(font);
        FontMetrics fm = g.getFontMetrics();
        for (Feat f : row) {
            int x1 = v.xOf(f.start - 1, LEFT);
            int x2 = v.xOf(f.end, LEFT);
            int bw = Math.max(2, x2 - x1);
            g.setColor(typeColor(f.type));
            g.fillRoundRect(x1, y, bw, h, 6, 6);
            if (f == hover) {
                g.setColor(Color.BLACK);
                g.drawRoundRect(x1, y, bw, h, 6, 6);
            }
            if (label && bw > 24) {
                String s = f.name;
                if (fm.stringWidth(s) + 8 < bw) {
                    g.setColor(Color.WHITE);
                    g.drawString(s, x1 + 6, y + h / 2 + fm.getAscent() / 2 - 2);
                }
            }
            if (bw > 10) {
                g.setColor(new Color(255, 255, 255, 180));
                int mid = y + h / 2;
                if (f.strand == '+') {
                    g.drawLine(x2 - 7, mid - 4, x2 - 2, mid);
                    g.drawLine(x2 - 7, mid + 4, x2 - 2, mid);
                } else if (f.strand == '-') {
                    g.drawLine(x1 + 7, mid - 4, x1 + 2, mid);
                    g.drawLine(x1 + 7, mid + 4, x1 + 2, mid);
                }
            }
        }
    }

    void label(Graphics2D g, String s, int y) {
        g.setFont(new Font("SansSerif", Font.BOLD, 11));
        g.setColor(new Color(90, 90, 95));
        g.drawString(s, 8, y);
    }

    static Color baseColor(char b) {
        switch (b) {
            case 'A': return browser.A;
            case 'T':
            case 'U': return browser.T;
            case 'G': return browser.G;
            case 'C': return browser.C;
            default: return browser.N;
        }
    }

    static Color typeColor(String t) {
        if (t.equals("gene")) return browser.GENE;
        if (t.equals("CDS")) return browser.CDS;
        if (t.contains("UTR")) return browser.UTR;
        return browser.OTHER;
    }

    static int tick(int span) {
        if (span <= 0) return 1;
        double raw = span / 8.0;
        double exp = Math.pow(10, Math.floor(Math.log10(raw)));
        double m = raw / exp;
        if (m < 2) return (int) exp;
        if (m < 5) return (int) (2 * exp);
        return (int) (5 * exp);
    }

    static String fmt(int p) {
        if (p >= 1000) return String.format("%,d", p);
        return Integer.toString(p);
    }

    Feat hit(int mx, int my) {
        int y0 = AXIS + SEQ + GAP;
        int a = v.visStart();
        int b = v.visEnd(trackW());
        for (Feat f : feats) {
            if (!f.overlaps(a, b)) continue;
            int lane = f.type.equals("gene") ? 0 : f.type.equals("CDS") ? 1 : 2;
            int y = y0 + lane * LANE;
            int h = f.type.equals("gene") ? 22 : f.type.equals("CDS") ? 16 : 12;
            int x1 = v.xOf(f.start - 1, LEFT);
            int x2 = v.xOf(f.end, LEFT);
            if (mx >= x1 && mx <= x2 && my >= y && my <= y + h) return f;
        }
        return null;
    }

    void updateStatus(int mx, int my) {
        if (status == null) return;
        int pos = (int) Math.floor(v.posOf(mx, LEFT)) + 1;
        if (pos < 1 || pos > fa.len() || mx < LEFT) {
            hover = null;
            status.setText(" drag to pan   wheel to zoom   Fit resets");
            repaint();
            return;
        }
        String s = " " + pos + " / " + fa.len() + "   " + fa.at(pos - 1);
        hover = hit(mx, my);
        if (hover != null) {
            s += "    " + hover.type + "  " + hover.name + "  "
                    + hover.start + "-" + hover.end + "  " + hover.strand;
        }
        status.setText(s);
        repaint();
    }

    public void mousePressed(MouseEvent e) {
        lastX = e.getX();
        setCursor(Cursor.getPredefinedCursor(Cursor.MOVE_CURSOR));
    }

    public void mouseReleased(MouseEvent e) {
        setCursor(Cursor.getDefaultCursor());
    }

    public void mouseDragged(MouseEvent e) {
        v.pan(e.getX() - lastX, trackW());
        lastX = e.getX();
        updateStatus(e.getX(), e.getY());
        repaint();
    }

    public void mouseMoved(MouseEvent e) {
        updateStatus(e.getX(), e.getY());
    }

    public void mouseWheelMoved(MouseWheelEvent e) {
        v.zoom(e.getX() - LEFT, e.getPreciseWheelRotation() > 0 ? 0.85 : 1.18);
        v.clamp(trackW());
        updateStatus(e.getX(), e.getY());
        repaint();
    }

    public void mouseClicked(MouseEvent e) {}
    public void mouseEntered(MouseEvent e) {}
    public void mouseExited(MouseEvent e) {
        hover = null;
        if (status != null) status.setText(" drag to pan   wheel to zoom   Fit resets");
        repaint();
    }
}

class Io {
    static BufferedReader open(String path) {
        try {
            InputStream in = new FileInputStream(path);
            if (path.endsWith(".gz")) in = new GZIPInputStream(in);
            return new BufferedReader(new InputStreamReader(in));
        } catch (IOException e) {
            browser.die("opening " + path, e);
            return null;
        }
    }
}
