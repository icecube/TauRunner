package web;
import java.awt.TextArea;
import java.io.OutputStream;
import java.awt.*;
import java.awt.geom.*;

/**
 * This class creates OutputStream object for TextArea, necessary for output with println.
 */

public class textStream extends OutputStream{
    private TextArea out=null;
    private byte a[];

    public textStream(TextArea out){
	this.out=out;
	a = new byte[1];
    }

    public void write(int b){
	a[0]=(byte) b;
	String str = new String(a);
	try{
	    if(out!=null) out.append(str);
	    if(aut!=null) gPrint(str);
	} catch(Exception e) {
	}
    }

    public void write(byte[] b){
	String str = new String(b);
	try{
	    if(out!=null) out.append(str);
	    if(aut!=null) gPrint(str);
	} catch(Exception e) {
	}
    }

    public void write(byte[] b, int off, int len){
	String str = new String(b, off, len);
	try{
	    if(out!=null) out.append(str);
	    if(aut!=null) gPrint(str);
	} catch(Exception e) {
	}
    }
    // almc begin

    private String[] buf;
    private String cur;
    private int pos=0;

    private Frame aut=null;
    private Rectangle2D.Double r2d=null;

    public textStream(Frame aut){
	this.aut=aut;
	r2d = new Rectangle2D.Double(100, 230, 400, 180);
	buf = new String[10];
	for(int i=0; i<10; i++) buf[i]="";
	a = new byte[1];
    }

    private synchronized void gPrint(String str){
	String sub;
	int end, len;

	sub=str;
	while(true){
	    len=sub.length(); if(len==0) break;
	    end=sub.indexOf("\n"); if(end==-1) end=len-1;
	    if(buf[pos].endsWith("\n")){ pos++; if(pos>=10) pos=0; buf[pos]=sub.substring(0, end+1); }
	    else buf[pos]+=sub.substring(0, end+1);
	    sub=sub.substring(end+1, len);
	}

	Graphics2D g2d = (Graphics2D)aut.getGraphics();
	g2d.setColor(Color.black);
	g2d.fill(r2d);
	g2d.setClip(r2d);
	g2d.setColor(Color.gray);
	for(int i=0; i<10; i++){
	    sub=buf[(pos+i+1)%10];
	    end=sub.indexOf("\n");
	    if(end!=-1) sub=sub.substring(0, end);
	    g2d.drawString(sub, (int)r2d.x, (int)(r2d.y+15+18*i));
	}
    }

    public void reflush(){
	gPrint("");
    }

    // almc end
}
