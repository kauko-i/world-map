import java.util.Date;
import java.util.Stack;
import java.util.ArrayList;
import java.util.Arrays;
import java.text.SimpleDateFormat;
import java.sql.*;
import java.io.File;
import javax.imageio.ImageIO;
import java.io.IOException;
import java.awt.*;
import java.awt.image.BufferedImage;
import java.awt.geom.AffineTransform;
import java.lang.Math;
import java.util.Scanner;
import java.util.Collections;

public class Kartta {

    public static final double MAAN_SADE = 6371;
    public static final double TARPEELLINEN_SADANTA = 1000;

    public static class Ruutu {
        private int i, j, kaytto, valtio, kasvukausi;
        private double lat, lon, ele, A, lammityskulu, rinteisyys, sadanta, kaikkiSade, keskilampo, jokivirtaa;
        private boolean jarvea, asuttu;
        public Ruutu(double lat, double lon, double ele, int i, int j) {
            this.lat = lat;
            this.lon = lon;
            this.ele = ele;
            this.i = i;
            this.j = j;
            this.A = Kartta.distance(lat - 0.25, lon, lat + 0.25, lon)*Kartta.distance(lat, lon - 0.25, lat, lon + 0.25);
        }

        public static double distance(Ruutu ruutu1, Ruutu ruutu2) {
            return Kartta.distance(ruutu1.lat, ruutu1.lon, ruutu2.lat, ruutu2.lon);
        }
    }

    public static class Saaasema {
        private double lat, lon;
        private double[] lampotila;
        private double[] sadanta;
        public Saaasema(double lat, double lon, double[] lampotila, double[] sadanta) {
            this.lat = lat;
            this.lon = lon;
            this.lampotila = lampotila;
            this.sadanta = sadanta;
        }
    }

    public static class Joki {
        private String nimi;
        private double virtaama;
        private ArrayList<Ruutu> reitti;
        public Joki(String nimi, double virtaama, ArrayList<Ruutu> reitti) {
            this.nimi = nimi;
            this.virtaama = virtaama;
            this.reitti = reitti;
        }
    }

	public static interface FunktioRuutuBoolean {
		public boolean f(Ruutu r);
	}
    
	public static interface FunktioRuutuVoid {
		public void f(Ruutu r);
	}

    public Ruutu[][] sisalto;

    public static double distance(double lat1, double lon1, double lat2, double lon2) {
        double lat1rad = lat1*Math.PI/180;
        double lon1rad = lon1*Math.PI/180;
        double lat2rad = lat2*Math.PI/180;
        double lon2rad = lon2*Math.PI/180;
        double latSin = Math.sin((lat2rad - lat1rad)/2);
        double lonSin = Math.sin((lon2rad - lon1rad)/2);
        double a = latSin*latSin + Math.cos(lat1rad)*Math.cos(lat2rad)*lonSin*lonSin;
        return MAAN_SADE*2*Math.atan2(Math.sqrt(a), Math.sqrt(1 - a));
    }

    public void piirra(int ruutuleveys) {
        BufferedImage bImg = new BufferedImage(ruutuleveys*this.sisalto.length, ruutuleveys*this.sisalto[0].length, BufferedImage.TYPE_INT_RGB);
        Graphics g = (Graphics2D)bImg.getGraphics();
        for (int i = 0; i < this.sisalto.length; i++) {
            for (int j = 0; j < this.sisalto[i].length; j++) {
                if (this.sisalto[i][j].ele == 0 || this.sisalto[i][j].jarvea) g.setColor(Color.blue);
                else if (this.sisalto[i][j].kaytto == 2) g.setColor(Color.orange);
                else if (this.sisalto[i][j].kaytto == 3) g.setColor(Color.green);
                else g.setColor(new Color(0,128,0));
                g.fillRect(i*ruutuleveys, (this.sisalto[i].length-j+1)*ruutuleveys, ruutuleveys, ruutuleveys);
                if (this.sisalto[i][j].asuttu) {
                    g.setColor(Color.red);
                    g.fillRect(i*ruutuleveys, (this.sisalto[i].length-j+1)*ruutuleveys, 1, ruutuleveys);
                }
            }
        }
        g.setColor(Color.white);
        Color[] varit = new Color[]{Color.blue, new Color(0,128,0), Color.green, Color.orange};
        String[] nimet = new String[]{"Vesistö", "Luonnonsuojelualue", "Talousmetsä", "Pelto", "Asuttu alue (väestötiheys ~1000/km2)"};
        g.fillRect(0, ruutuleveys*this.sisalto[0].length - 15*(varit.length + 1) - 5, 300, ruutuleveys*this.sisalto[0].length);
        for (int i = 0; i < varit.length + 1; i++) {
            if (i < varit.length) {
                g.setColor(varit[i]);
                g.fillRect(10, ruutuleveys*this.sisalto[0].length - 15*(varit.length + 1) + i*15, 10, 10);
            }
            g.setColor(Color.black);
            g.drawString(nimet[i], 25, ruutuleveys*this.sisalto[0].length - 15*(varit.length + 1) + i*15 + 10);
        }
        g.setColor(Color.red);
        for (int i = 10; i < 20; i += 2) g.fillRect(i, ruutuleveys*this.sisalto[0].length - 15, 1, 10);
        try {
            if (ImageIO.write(bImg, "png", new File("kartat/"+new SimpleDateFormat("yyMMddHHmm").format(new Date())+".png"))) {
                System.out.println("Saved");
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public boolean kartalla(int i, int j) {
        return -1 < i && -1 < j && i < this.sisalto.length && j < this.sisalto[i].length;
    }

    public void floodfill(int i, int j, FunktioRuutuBoolean ehto, FunktioRuutuVoid mita) {
        Stack<Ruutu> pino = new Stack<Ruutu>();
        pino.push(this.sisalto[i][j]);
        boolean[][] vierailtu = new boolean[this.sisalto.length][this.sisalto[0].length];
        while (!pino.empty()) {
            Ruutu next = pino.pop();
            if (vierailtu[next.i][next.j]) continue;
            vierailtu[next.i][next.j] = true;
            if (!ehto.f(next)) continue;
            mita.f(next);
            if (kartalla(next.i-1, next.j) && !vierailtu[next.i-1][next.j]) pino.push(this.sisalto[next.i-1][next.j]);
            if (kartalla(next.i+1, next.j) && !vierailtu[next.i+1][next.j]) pino.push(this.sisalto[next.i+1][next.j]); 
            if (kartalla(next.i, next.j-1) && !vierailtu[next.i][next.j-1]) pino.push(this.sisalto[next.i][next.j-1]); 
            if (kartalla(next.i, next.j+1) && !vierailtu[next.i][next.j+1]) pino.push(this.sisalto[next.i][next.j+1]); 
        }
    }

    public void bresenham(int x1, int y1, int x2, int y2, FunktioRuutuVoid kasittely) {
        int minx = Math.min(x1, x2);
        int miny = Math.min(y1, y2);
        int maxx = Math.max(x1, x2);
        int maxy = Math.max(y1, y2);
        double d = 1.0*(y2 - y1)/(x2 - x1);
        if (maxy - miny < maxx - minx) {
            for (int i = minx; i <= maxx; i++) {
                int j = (int)((i - x1)*d + y1);
                if (kartalla(i, j)) kasittely.f(sisalto[i][j]);
            }
        } else {
            for (int j = miny; j <= maxy; j++) {
            	int i = (int)((j - y1)/d + x1);
                if (kartalla(i, j)) kasittely.f(sisalto[i][j]);
            }
        }
    }

    public static void main(String[] args) throws Exception {
        Kartta map = new Kartta();
        map.sisalto = new Ruutu[720][357];
        System.out.println("Hello World!!!");
        Class.forName("org.sqlite.JDBC");
        Connection conn = DriverManager.getConnection("jdbc:sqlite:karttadata.db");
        Statement stat = conn.createStatement();

        ResultSet rs = stat.executeQuery("select * from ruudut;");
        while (rs.next()) {
            double lat = Double.parseDouble(rs.getString("lat"));
            double lon = Double.parseDouble(rs.getString("lon"));
            double ele = Double.parseDouble(rs.getString("ele"));
            int i = (int)((lon + 180)*2);
            int j = (int)((lat + 89)*2);
            Ruutu uusiRuutu = new Ruutu(lat, lon, ele, i, j);
            map.sisalto[i][j] = uusiRuutu;
        }
        stat.executeQuery("select * from joki;");
        ArrayList<Joki> joet = new ArrayList<Joki>();
        while (rs.next()) {
            String joki = rs.getString("nimi");
            double virtaama = rs.getDouble("virtaama");
            joet.add(new Joki(joki, virtaama, new ArrayList<Ruutu>()));
        }
        for (Joki joki : joet) {
            stat.executeQuery(String.format("select * from jokimutka where joki=\"%s\" order by nro", joki.nimi));
            Ruutu edellinen = null;
            while (rs.next()) {
                double lat = rs.getDouble("lat");
                double lon = rs.getDouble("lon");
                int i = (int)((lon + 180)*2);
                int j = (int)((lat + 89)*2);
                if (edellinen != null) {
                    map.bresenham(edellinen.i, edellinen.j, i, j, r -> r.jokivirtaa += joki.virtaama);
                }
                edellinen = map.sisalto[i][j];
            }
        }
        rs.close();
        conn.close();
        for (int i = 0; i < map.sisalto.length; i++) {
            for (int j = 0; j < map.sisalto[i].length; j++) {
                if (map.sisalto[i][j].ele != 0) {
                    double ele = map.sisalto[i][j].ele;
                    ArrayList<Ruutu> jarvi = new ArrayList<Ruutu>();
                    map.floodfill(i, j, r -> r.ele == ele, r -> jarvi.add(r));
                    if (4 < jarvi.size()) {
                        for (Ruutu r : jarvi) r.jarvea = true;
                    }
                }
            }
        }
        Scanner ilmastolukija = new Scanner(new File("saahistoria.dat"));
        ArrayList<Saaasema> saaasemat = new ArrayList<Saaasema>();
        while (ilmastolukija.hasNextLine()) {
            String[] asemadata = ilmastolukija.nextLine().split("\\|");
            String[] lammot = asemadata[2].split(",");
            String[] sadannat = asemadata[3].split(",");
            double[] lammotDouble = new double[lammot.length];
            double[] sadannatDouble = new double[sadannat.length];
            for (int i = 0; i < lammot.length; i++) {
                lammotDouble[i] = Double.parseDouble(lammot[i]);
                sadannatDouble[i] = Double.parseDouble(sadannat[i]);
            }
            saaasemat.add(new Saaasema(Double.parseDouble(asemadata[0]), Double.parseDouble(asemadata[1]), lammotDouble, sadannatDouble));
        }
        ilmastolukija.close();
        double maxLammityskulu = 0;
        double maxsade = 0;
        double maxsadanta = 0;
        double maxlampo = Double.NEGATIVE_INFINITY;
        double minlampo = Double.POSITIVE_INFINITY;
        for (int i = 0; i < map.sisalto.length; i++) {
            System.out.println(i);
            for (int j = 0; j < map.sisalto[i].length; j++) {
                final double LAT = map.sisalto[i][j].lat;
                final double LON = map.sisalto[i][j].lon;
                ArrayList<Saaasema> ehdokasasemat = new ArrayList<Saaasema>(saaasemat);
                ehdokasasemat.removeIf(a -> 5 < Math.abs(a.lat - LAT) || 5 < Math.abs(a.lon - LON));
                Collections.sort(ehdokasasemat, (a1, a2) -> {
                    double e1 = distance(LAT, LON, a1.lat, a1.lon);
                    double e2 = distance(LAT, LON, a2.lat, a2.lon);
                    if (e1 == e2) return 0;
                    return e1 < e2 ? -1 : 1;
                });
                double[] lampoarvio = new double[365];
                double[] sadearvio = new double[365];
                double sadanta = 0;
                double lammityskulu = 0;
                double kaikkiSade = 0;
                double keskilampo = 0;
                int kasvukausi = 0;
                for (int k = 0; k < lampoarvio.length; k++) {
                    if (distance(map.sisalto[i][j].lat, map.sisalto[i][j].lon, ehdokasasemat.get(0).lat, ehdokasasemat.get(0).lon) == 0) {
                        lampoarvio[k] = ehdokasasemat.get(0).lampotila[k];
                        sadearvio[k] = ehdokasasemat.get(0).sadanta[k];
                    } else {
                        double dsumma = 0;
                        for (int l = 0; l < Math.min(ehdokasasemat.size(), 3); l++) {
                            double d = distance(map.sisalto[i][j].lat, map.sisalto[i][j].lon, ehdokasasemat.get(l).lat, ehdokasasemat.get(l).lon);
                            lampoarvio[k] += ehdokasasemat.get(l).lampotila[k]/d;
                            sadearvio[k] += ehdokasasemat.get(l).sadanta[k]/d;
                            dsumma += 1/d;
                        }
                        lampoarvio[k] /= dsumma;
                        sadearvio[k] /= dsumma;
                    }
                    lammityskulu += Math.max(0, 15 - lampoarvio[k]) + 2*Math.max(0, lampoarvio[k] - 30);
                    sadanta += Math.max(3, sadearvio[k])*Math.pow(Math.max(100 - lampoarvio[k], 0)/100, 3);
                    kaikkiSade += sadearvio[k]*Math.pow(Math.max(100 - lampoarvio[k], 0)/100, 3);
                    keskilampo += lampoarvio[k];
                    if (10 < lampoarvio[k]) kasvukausi++;
                }
                map.sisalto[i][j].lammityskulu = lammityskulu;
                maxLammityskulu = Math.max(maxLammityskulu, lammityskulu);
                map.sisalto[i][j].sadanta = Math.min(sadanta, TARPEELLINEN_SADANTA);
                map.sisalto[i][j].kaikkiSade = kaikkiSade;
                maxsade = Math.max(maxsade, kaikkiSade);
                map.sisalto[i][j].keskilampo = keskilampo/lampoarvio.length;
                maxlampo = Math.max(maxlampo, keskilampo/lampoarvio.length);
                minlampo = Math.min(minlampo, keskilampo/lampoarvio.length);
                map.sisalto[i][j].kasvukausi = kasvukausi;
            }
        }
        double maxRinteisyys = 0;
        ArrayList<Ruutu> maa = new ArrayList<Ruutu>();
        double[][] rinteisyys = new double[map.sisalto.length][map.sisalto[0].length];
        for (int i = 0; i < map.sisalto.length; i++) {
            System.out.println(i);
            for (int j = 0; j < map.sisalto[i].length; j++) {
                if (map.sisalto[i][j].ele != 0 && !map.sisalto[i][j].jarvea) {
                    maa.add(map.sisalto[i][j]);
                    ArrayList<int[]> suunnat = new ArrayList<int[]>();
                    suunnat.add(new int[]{-1,0});
                    suunnat.add(new int[]{1,0});
                    suunnat.add(new int[]{0,-1});
                    suunnat.add(new int[]{0,1});
                    int naapureita = 0;
                    for (int[] suunta : suunnat) {
                        if (map.kartalla(i + suunta[0], j + suunta[1])) {
                            Ruutu naapuri = map.sisalto[i + suunta[0]][j + suunta[1]];
                            if (naapuri.ele != 0 && !naapuri.jarvea) {
                                naapureita++;
                                rinteisyys[i][j] += (naapuri.ele - map.sisalto[i][j].ele)*(naapuri.ele - map.sisalto[i][j].ele);
                            }
                        }
                    }
                    if (naapureita != 0) {
                        map.sisalto[i][j].rinteisyys = rinteisyys[i][j]/naapureita;
                        maxRinteisyys = Math.max(maxRinteisyys, rinteisyys[i][j]/naapureita);
                    } else map.sisalto[i][j].rinteisyys = 0;
                }
            }
        }
        final double MAX_SADE = maxsade;
        final double MAX_LAMPO = maxlampo;
        final double MIN_LAMPO = minlampo;
        Collections.sort(maa, (r1, r2) -> {
            double r1x = 2*r1.kaikkiSade/MAX_SADE + (r1.keskilampo - MIN_LAMPO)/(MAX_LAMPO - MIN_LAMPO);
            double r2x = 2*r2.kaikkiSade/MAX_SADE + (r2.keskilampo - MIN_LAMPO)/(MAX_LAMPO - MIN_LAMPO);
            if (r1x == r2x) return 0;
            return r1x < r2x ? 1 : -1;
        });
        double tropiikkia = 0;
        int i = 0;
        while (tropiikkia < 10000000) {
            Ruutu next = maa.get(i++);
            next.kaytto = 1;
            tropiikkia += next.A;
        }
        maa = new ArrayList<Ruutu>(maa.subList(i, maa.size()));
        final double MAX_LAMMITYSKULU = maxLammityskulu;
        final double MAX_RINTEISYYS = maxRinteisyys;
        Collections.sort(maa, (r1, r2) -> {
            double r1x = 2*((r1.sadanta + r1.jokivirtaa)/TARPEELLINEN_SADANTA) - r1.lammityskulu/MAX_LAMMITYSKULU - Math.sqrt(r1.rinteisyys/MAX_RINTEISYYS);
            double r2x = 2*((r2.sadanta + r2.jokivirtaa)/TARPEELLINEN_SADANTA) - r2.lammityskulu/MAX_LAMMITYSKULU - Math.sqrt(r2.rinteisyys/MAX_RINTEISYYS);
            if (r1x == r2x) return 0;
            return r1x < r2x ? 1 : -1;
        });
        double asuttua = 0;
        i = 0;
        while (asuttua < 12000000) {
            Ruutu next = maa.get(i++);
            next.asuttu = true;
            asuttua += next.A;
        }
        Collections.sort(maa, (r1, r2) -> {
            double r1x = 2*((r1.sadanta + r1.jokivirtaa)/TARPEELLINEN_SADANTA) + Math.min(180, r1.kasvukausi)/180 - 2*Math.sqrt(r1.rinteisyys/MAX_RINTEISYYS);
            double r2x = 2*((r2.sadanta + r2.jokivirtaa)/TARPEELLINEN_SADANTA) + Math.min(180, r2.kasvukausi)/180 - 2*Math.sqrt(r2.rinteisyys/MAX_RINTEISYYS);
            if (r1x == r2x) return 0;
            return r1x < r2x ? 1 : -1;
        });
        double peltoa = 0;
        i = 0;
        while (peltoa < 12000000) {
            Ruutu next = maa.get(i++);
            next.kaytto = 2;
            peltoa += next.A;
        }
        Collections.sort(maa, (r1, r2) -> {
            double r1x = 2*r1.kaikkiSade/MAX_SADE + (r1.keskilampo - MIN_LAMPO)/(MAX_LAMPO - MIN_LAMPO);
            double r2x = 2*r2.kaikkiSade/MAX_SADE + (r2.keskilampo - MIN_LAMPO)/(MAX_LAMPO - MIN_LAMPO);
            if (r1x == r2x) return 0;
            return r1x < r2x ? 1 : -1;
        });
        double talousmetsaa = 0;
        i = 0;
        while (talousmetsaa < 3000000) {
            Ruutu next = maa.get(i++);
            next.kaytto = 3;
            talousmetsaa += next.A;
        }
        map.piirra(2);
    }
}
